##### Performance of LSMM when the random effects are correlated #####
# Vary rho=0, 0.2, 0.4, 0.6, 0.8 to get Figure S52-53 in Supplementary Document

library(LSMM)
library(pROC)
library(MASS)

source("performance.R")

# function to generate data
generate_data_corr_beta <- function(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, rho){
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
  Sigma <- sigma2*rho^(abs(outer(1:K, 1:K, "-")))
  beta <- mvrnorm(1, rep(0, K), Sigma)
  beta[-indexeta] <- 0
  
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
K      <- 500      # No. of random effects
Z.perc <- 0.1      # the proportion the entries in Z is 1
A.perc <- 0.1      # the proportion the entries in A is 1
alpha  <- 0.2      # parameter in the Beta distribution
beta0  <- -2       # intercept of the logistic model
set.seed(1)
b      <- rnorm(L) # fixed effects
omega  <- 0.2      # proportion of relevant annotations
sigma2 <- 1        # parameter in the spike-slab prior
rho    <- 0.2      # rho in the autocorrelation matrix
rep    <- 50       # repeat times

result <- matrix(0, rep, 16)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")
  
  data <- generate_data_corr_beta(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, corr)
  
  fit <- LSMM(data$Pvalue, data$Z, data$A)
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")
  
  result[i, 1:4]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma, 1-fit$pi1))
  result[i, 5:8]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma.stage1, 1-fit$pi1.stage1))
  result[i, 9:12] <- as.numeric(performance(data$gamma, assoc.SNP$gamma.stage2, 1-fit$pi1.stage2))

  relev.Anno <- relev.Anno(fit, FDRset = 0.1, fdrControl = "global")
  result[i, 13:16] <- as.numeric(performance(data$eta, relev.Anno, 1-fit$omegak))
}

result <- as.data.frame(result)
names(result) <- c("FDR.LSMM", "power.LSMM", "AUC.LSMM", "pAUC.LSMM", 
                   "FDR.TGM", "power.TGM", "AUC.TGM", "pAUC.TGM", 
                   "FDR.LFM", "power.LFM", "AUC.LFM", "pAUC.LFM",
                   "FDR.Anno", "power.Anno", "AUC.Anno", "pAUC.Anno")
