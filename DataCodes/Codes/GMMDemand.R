########################################
# Código para rodar o modelo de demanda
# Escolha Discreta com Dados Agregados
# Claudio R. Lucinda
# USP
# 2020
########################################
# SUGEST?O: Atualizar pra 3.6.0 o R
########################################

#install.packages("reticulate")
#install.packages("tidyverse")
#install.packages("SASxport")
#install.packages("Hmisc")
#install.packages("AER")
#install.packages("ivreg")
#install.packages("gmm")
#install.packages("plot.matrix")


library("reticulate")
library("tidyverse")
library("Hmisc")
library("SASxport")
library("AER")
library("gmm")
library("plot.matrix")

dat<-sasxport.get("./Data/NevoData_OI.xpt")

dat<-dat %>%
  rename(price=v2, share=v1, sugar=v29, mushy=v30, cdid=v75, id=v76) 

data<-dat[,c("price", "share", "sugar", "mushy", "cdid", "id")]

data<-data %>%
  group_by(cdid) %>%
  mutate(inshare=sum(share))

data$outshare<-1-data$inshare

data$meanu<-log(data$share)-log(data$outshare)


# Basicão - OLS
model_0<-lm(meanu~price+sugar+mushy, data=data)
summary(model_0)

# Basicão - OLS + Efeitos Fixos de região
model_1<-lm(meanu~price+sugar+mushy+factor(cdid), data=data)
summary(model_1)

# Criando as dummies de produto
dum_data<-dat[,4:27]
teste <- simplify2array(
  apply(
    dum_data, 1, 
    function(x) paste(names(dum_data)[x != 0], collapse = " ")
  )
)
testevar<-data.frame(teste)
testevar$teste<-as.factor(testevar$teste)

data02<-data
data02$brand_id<-testevar$teste

# Dummies de Produto
model_2<-lm(meanu~price+sugar+mushy+factor(cdid)+factor(brand_id), data=data02)
summary(model_2)


# IVReg
instruments<-data.frame(dat[,32:75])
data03<-data.frame(data02,instruments)
inst_names<-colnames(instruments)
fla<-paste0("~",paste(inst_names, collapse="+"),"+sugar+mushy")

model_3<-ivreg(meanu~price+sugar+mushy,instruments=as.formula(fla), data=data03)
summary(model_3, diagnostics = TRUE)

set.seed(1024)
model_3_gmm <- gmm(meanu~price+sugar+mushy, # linear regression formula
                   as.formula(fla), # formular to specify the instruments and exogenous variables
                   rnorm(length(coef(model_3))),
                   type="twoStep",
                   data = data03,
                   wmatrix = "ident",
                   vcov = "MDS",
                   optfct = "optim",
                   control = list(eval.max = 10000))
summary(model_3_gmm)
