##############################################
# Estimação GMM e ML
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
# Poisson (variance equals mean) with Default standard errors based on the
# inverse of the Hessian of log like function

poisson_res <- patents_df %>%
  glm(p91 ~ lr91 + aerosp + chemist + computer + machines +
        vehicles + japan + us, family = poisson(), data = .)

summary(poisson_res)

# GMM estimator as Quasi- (Psudo-) MLE
# Define moment conditions matrix for QMLE (gmm)

mom <- function(beta, df) {
  # df is the data frame with first column as dv
  # This function returns n * q matrix
  # Each column is one moment condition before taking sample average
  # There are totally q moment conditions
  
  y <- as.numeric(df[, 1])
  x <- data.matrix(df[, 2:ncol(df)])
  
  # Refer to moment conditions of QMLE
  m <- x * as.vector(y - exp(x %*% beta))
  
  return(cbind(m))
}

# Generate regression coef as the initial values for QMLE (gmm)
init_values <- patents_df %>%
  select(
    p91, lr91, aerosp, chemist, computer, machines,
    vehicles, japan, us
  ) %>%
  mutate(const = 1) %>%
  # Please hold the order as previous glm() to facilitate comparison
  select(
    p91, const, lr91, aerosp, chemist, computer, machines,
    vehicles, japan, us
  ) %>%
  lm(p91 ~ lr91 + aerosp + chemist + computer + machines +
       vehicles + japan + us, data = .) |>
  {\(x) coef(x)*0.01}()

# Be careful that we need to use "nlminb" INSTEAD OF "optim", which is
# bloody awful.
# Refer to https://cran.r-project.org/web/packages/gmm/vignettes/gmm_with_R.pdf

# Just use random initial values
set.seed(1024)
qmle_res_1 <- gmm(mom, df_n, rnorm(length(init_values)),
                  wmatrix = "optimal",
                  vcov = "MDS",
                  optfct = "nlminb",
                  control = list(eval.max = 10000)
)

# Use initial values from regression coef
set.seed(1024)
qmle_res_2 <- gmm(mom, df_n, init_values,
                  wmatrix = "optimal",
                  vcov = "MDS",
                  optfct = "nlminb",
                  control = list(eval.max = 10000)
)

summary(qmle_res_1)
summary(qmle_res_2)