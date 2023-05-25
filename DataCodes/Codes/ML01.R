#########################################################
# Código para Replicar alguns dos exemplos do Greene
# Claudio R. Lucinda
# 2023
# FEA-USP
##########################################################

set.seed(2512)

library(tidyverse)
library(ggplot2)

# Pegando os dados do Exemplo C.1 do Greene

Data=read.csv("./Data/TableFC-1.csv")

# Colocando a função log-verossimilhança
# Passo 1 - Calcular na mão o vetor de contribuições de log-verossimilhança
# Dá uma ideia boa de coisas que podem dar errado

beta=.1

loglik=-(log(beta+Data$E)+(Data$Y/(beta+Data$E)))
sum(loglik)

# Passo 2 - Fazer a verossimilhança como uma função
# Prestar atenção à diferença de nomes aqui. 

llik <- function(parm, data) {
  
  loglik=-(log(parm+data$E)+(Data$Y/(parm+data$E)))
  return(sum(loglik))
  
}



llik(parm=.1, data=Data)


# Deu igual, portanto beleza

# Vamos fazer um gráfico jeitosinho da função - como o vetor de parâmetros 
# Lógico que a forma e a posição disso vai depender dos dados

points=seq(from=1,to=30, length.out=200)

logliks = map_dbl(points, llik, data=Data)

plot(points, logliks,main = "-Log Likelihood",
     xlab = "Parametro", ylab = "Log Likelihood",
     pch = 19, frame = FALSE)


# Vamos achar o estimador ML do beta
# Bem como outras coisinhas
# Como o optim vai procurar o mínimo, vou criar a função "virada"


llik2 <- function(parm, data) {
  
  loglik=(log(parm+data$E)+(Data$Y/(parm+data$E)))
  return(sum(loglik))
  
}

# Vamos minimizar isso
# Vários estimadores - vai dar tudo igual pois LL "não cabulosa"
start = 10
# ThetaNM = optim(start,llik2)
ThetaSANN = optim(start,llik2, data = Data, method="SANN")
ThetaBFGS = optim(start,llik2, data = Data, hessian=T,method="BFGS")

SeBFGS = sqrt(diag(solve(ThetaBFGS$hess)))

# Fazendo os testes de hipóteses
# O que temos acima é a versão restrita de um caso mais geral
# ver p. 490 Greene


llik3 <- function(parm, data) {
  beta=1/(parm[1]+data$E)
  loglik=-log(((beta^parm[2])/gamma(parm[2]))*(data$Y^(parm[2]-1))*exp(-data$Y*beta))
  return(sum(loglik))
}

start2 = c(-3,5)

ThetaBFGS2 = optim(start2,llik3, data = Data, hessian=T,method="BFGS")

SeBFGS2 = sqrt(diag(solve(ThetaBFGS2$hess)))

# Testando a restrição (rho=1)
# LR
LRStat=-2*(-ThetaBFGS$value+ThetaBFGS2$value)

#Wald
WaldStat=((ThetaBFGS2$par[2]-1)^2)/(SeBFGS2[2]^2)

# LM
# Preciso dos gradientes

grad <- function(parm, data) {
  beta=1/(parm[1]+data$E)
  n=length(data$E)
  dldbeta=sum(-parm[2]*beta+data$Y*beta^2)
  dldrho=sum(log(beta))-n*digamma(parm[2])+sum(log(data$Y))
  return(c(dldbeta, dldrho))
  
}

hess <- function(parm, data) {
  beta=1/(parm[1]+data$E)
  n=length(data$E)
  d2ldbeta2=sum(parm[2]*(beta^2)-2*data$Y*beta^3)
  d2ldrho2=-n*trigamma(parm[2])
  d2ldbetadrho=-sum(beta)
  return(matrix(c(d2ldbeta2,d2ldbetadrho, d2ldbetadrho, d2ldrho2), nrow=2, ncol=2))
}


GradRes=grad(c(ThetaBFGS$par[1],1), data=Data)
HessRes=hess(c(ThetaBFGS$par[1],1), data=Data)


LMStat=-t(GradRes)%*%solve(HessRes)%*%GradRes
