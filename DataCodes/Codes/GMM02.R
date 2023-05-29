##############################################
# Estimação GMM
# Claudio R. Lucinda
# 2023
# FEA/USP
##############################################

# Dados Simulados
# semente
#install.packages("gmm")
library(gmm)
library(tidyverse)
library(haven)

patents_df <- read_dta("./Data/patents.dta")
# explore
glimpse(patents_df)

# Fazendo Grafico

patents_df %>%
  ggplot(aes(p91)) %>%
  +geom_histogram()

# OLSzão da massa

lm_res <- patents_df %>%
  lm(p91 ~ lr91 + aerosp + chemist + computer + machines +
       vehicles + japan + us, data = .)

summary(lm_res)

# Ajeitando os dados

df_n <- patents_df %>%
  select(
    p91, lr91, aerosp, chemist, computer, machines,
    vehicles, japan, us
  ) %>%
  mutate(const = 1) %>%
  # Please hold the order as previous glm() to facilitate comparison
  select(
    p91, const, lr91, aerosp, chemist, computer, machines,
    vehicles, japan, us
  )

## Need to converse the tibble class to dataframe
df_n <- as.data.frame(df_n)

set.seed(1024)
gmm_lm_iid_1 <- gmm(
  p91 ~ lr91 + aerosp + chemist + computer + machines +
    vehicles + japan + us, # linear regression formula
  ~ lr91 + aerosp + chemist + computer + machines +
    vehicles + japan + us, # formular to specify the instruments and exogenous variables
  rnorm(length(coef(lm_res))),
  data = df_n,
  wmatrix = "optimal",
  vcov = "iid",
  optfct = "nlminb",
  control = list(eval.max = 10000)
)
summary(gmm_lm_iid_1)

#Twostep

set.seed(1024)
gmm_lm_iid_2 <- gmm(
  p91 ~ lr91 + aerosp + chemist + computer + machines +
    vehicles + japan + us, # linear regression formula
  df_n[, 2:length(df_n)], # matrix to specify the instruments and exogenous variables
  rnorm(length(coef(lm_res))),
  data = df_n,
  wmatrix = "optimal",
  vcov = "iid",
  optfct = "nlminb",
  control = list(eval.max = 10000)
)
summary(gmm_lm_iid_2)

#GMM Cov
diag(mean(resid(lm_res)^2) * solve(t(as.matrix(df_n[, 2:length(df_n)])) %*% as.matrix(df_n[, 2:length(df_n)])))^0.5

# LM Cov
diag(sum(resid(lm_res)^2)/(nrow(df_n)-length(df_n)+1) * solve(t(as.matrix(df_n[, 2:length(df_n)])) %*% as.matrix(df_n[, 2:length(df_n)])))^0.5


# Pra aproveitar toda a flexibilidade do comando, melhor usar a versão geral, ao invés da versão com fórmula

mom_lm <- function(beta, df) {
  # df is the data frame with first column as dv
  # This function returns n * q matrix
  # Each column is one moment condition before taking sample average
  # There are totally q moment conditions
  
  y <- as.numeric(df[, 1])
  x <- data.matrix(df[, 2:ncol(df)])
  
  # Refer to moment conditions of QMLE
  m <- x * as.vector(y - x %*% beta)
  
  return(cbind(m))
}


set.seed(1024)
gmm_lm_mds <- gmm(mom_lm, df_n, rnorm(length(coef(lm_res))),
                  wmatrix = "optimal",
                  vcov = "MDS",
                  optfct = "nlminb",
                  control = list(eval.max = 10000)
)
summary(gmm_lm_mds)

# Agora batem com a versão robusta do OLS
diag(vcovHC(lm_res, type = "HC0"))^0.5