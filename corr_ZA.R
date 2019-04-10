##### Performance in identification of relevant annotations when the annotations #####
##### in design matrix of fixed effects and random effects are not independent #####
# Vary corr=0, 0.2, 0.4, 0.6, 0.8 to get 
# Figure S23 in Supplementary Document

library(LSMM)
library(pROC)
library(MASS)

source("performance.R")

# function to generate data
generate_data_corr_ZA <- function(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, q){
  # design matrix of fixed effects & random effects
  corr <- matrix(0, K+L, K+L)
  corr[1:(L+50), 1:(L+50)] <- q
  diag(corr) <- 1
  preZA <- mvrnorm(M, rep(0, K+L), corr)
  ZA <- t(t(preZA) < apply(preZA, 2, quantile, probs = A.perc))
  Z <- ZA[, 1:L]
  A <- ZA[, (L+1):(K+L)]
  
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
  
  # Pvalue
  Pvalue             <- runif(M)
  Pvalue[indexgamma] <- rbeta(sum(indexgamma), alpha, 1)
  
  return( list(Z = Z, A = A, Pvalue = Pvalue, beta = beta, pi1 = pi1, eta = eta,
               gamma = gamma))
}

# sigmoid function
sigma <- function(x){
  y <- 1/(1+exp(-x))
  return (y)
}

M      <- 100000   # No. of SNPs
L      <- 10       # No. of fixed effects
K      <- 100      # No. of random effects
Z.perc <- 0.1      # the proportion the entries in Z is 1
A.perc <- 0.1      # the proportion the entries in A is 1
alpha  <- 0.2      # parameter in the Beta distribution
beta0  <- -2       # intercept of the logistic model
set.seed(1)
b      <- rnorm(L) # fixed effects
omega  <- 0.2      # proportion of relevant annotations
sigma2 <- 1        # parameter in the spike-slab prior
corr   <- 0.2      # correlation coefficient among the annotations with fixed effects and the first 50 the annotations with random effects
rep    <- 50       # repeat times

result <- matrix(0, rep, 8)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")
  
  data <- generate_data_corr_ZA(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, corr)

  # LSMM
  fit <- LSMM(data$Pvalue, data$Z, data$A)
  relev.Anno <- relev.Anno(fit, FDRset = 0.1, fdrControl = "global")
  result[i, 1:4] <- as.numeric(performance(data$eta, relev.Anno, 1-fit$omegak))
  
  # LSMM without fixed effects
  fit1 <- LSMM(data$Pvalue, NULL, data$A)
  relev.Anno1 <- relev.Anno(fit1, FDRset = 0.1, fdrControl = "global")
  result[i, 5:8] <- as.numeric(performance(data$eta, relev.Anno1, 1-fit$omegak))
}

result <- as.data.frame(result)
names(result) <- c("FDR.LSMM", "power.LSMM", "AUC.LSMM", "pAUC.LSMM",
                   "FDR.LSMM.no.Z", "power.LSMM.no.Z", "AUC.LSMM.no.Z", "pAUC.LSMM.no.Z")
