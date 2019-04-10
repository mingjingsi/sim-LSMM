##### Estimation of alpha using individual-level data#####
# Vary N=5000, 10000, 
# h2=0.2, 0.4, 0.6, 0.8 to get 
# Figures S40 in Supplementary Document

library(LSMM)
library(MASS)

# function to generate data
generate_data_ind <- function(M, N, L, K, Z.perc, A.perc, beta0, b, omega, sigma2, h2){
  
  # genotype data
  X    <- matrix(1, N, M)
  X1   <- matrix(runif(M*N), N, M)
  f    <- runif(M, 0.05, 0.5)
  f2   <- matrix(rep(f^2, each = N), N, M)
  f1_2 <- matrix(rep((1 - f)^2, each = N), N, M)
  X[X1 < f2]         <- 2
  X[(1 - X1) < f1_2] <- 0
  
  # design matrix of fixed effects
  Z         <- rep(0, M*L)
  indexZ    <- sample(M*L, M*L*Z.perc)
  Z[indexZ] <- 1
  Z         <- matrix(Z, M, L)
  
  # design matrix of random effects
  A         <- rep(0, M*K)
  indexA    <- sample(M*K, M*K*A.perc)
  A[indexA] <- 1
  A         <- matrix(A, M, K)
  
  # eta (latent variable which indicate whether the annotation is relevant to the phenotype)
  eta           <- rep(0, K)
  indexeta      <- sample(K, K*omega)
  eta[indexeta] <- 1
  
  # beta (random effects)
  beta           <- rep(0, K)
  beta[indexeta] <- rnorm(K*omega, 0, sqrt(sigma2))
  
  # gamma (latent variable which indicate whether the SNP is associated with the phenotype)
  pi1               <- sigma(beta0 + Z %*% b + A %*% beta)
  gamma             <- rep(0, M)
  indexgamma        <- (runif(M) < pi1)
  gamma[indexgamma] <- 1
  
  # effect size
  beta_SNP <- numeric(M)
  beta_SNP[indexgamma] <- rnorm(sum(indexgamma), 0, 1)
  
  # environment effect
  e <- rnorm(N, 0, sqrt((1/h2-1)*var(X%*%beta_SNP)))
  
  # phenotype
  y <- X%*%beta_SNP + e
  
  # p-value
  Pvalue <- t(apply(X, 2, function(X.col) summary(lm(y ~ X.col))$coefficients[2,4]))
  
  return( list(Pvalue = Pvalue, Z = Z, A = A, beta = beta, eta = eta, gamma = gamma))
}

# sigmoid function
sigma <- function(x){
  y <- 1/(1+exp(-x))
  return (y)
}

M      <- 20000    # No. of SNPs
N      <- 10000    # No. of individuals
L      <- 10       # No. of fixed effects
K      <- 100      # No. of random effects
Z.perc <- 0.1      # the proportion the entries in Z is 1
A.perc <- 0.1      # the proportion the entries in A is 1
beta0  <- -2       # intercept of the logistic model
set.seed(1)
b      <- rnorm(L) # fixed effects
omega  <- 0.1      # proportion of relevant annotations
sigma2 <- 1        # parameter in the spike-slab prior
h2     <- 0.8      # heritability
rep    <- 50       # repeat times

result <- numeric(rep)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")
  
  data <- generate_data_ind(M, N, L, K, Z.perc, A.perc, beta0, b, omega, sigma2, h2)
  
  fit <- LSMM(data$Pvalue, data$Z, data$A)
  
  result[i] <- fit$alpha
}

result <- as.data.frame(result)
names(result) <- c("est.alpha.ind")
