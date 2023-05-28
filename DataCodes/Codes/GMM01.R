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

set.seed(123)
# Dados Gerados da Distribuição normal
n <- 200
x <- rnorm(n, mean = 4, sd = 2)
# set up the moment conditions for comparison
# MM (just identified)
g0 <- function(tet, x) {
  m1 <- (tet[1] - x)
  m2 <- (tet[2]^2 - (x - tet[1])^2)
  f <- cbind(m1, m2)
  return(f)
}
# GMM (over identified)
g1 <- function(tet, x) {
  m1 <- (tet[1] - x)
  m2 <- (tet[2]^2 - (x - tet[1])^2)
  m3 <- x^3 - tet[1] * (tet[1]^2 + 3 * tet[2]^2)
  f <- cbind(m1, m2, m3)
  return(f)
}
print(res0 <- gmm(g0, x, c(mu = 0, sig = 0)))
print(res1 <- gmm(g1, x, c(mu = 0, sig = 0)))

summary(res0)