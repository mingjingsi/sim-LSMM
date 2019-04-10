##### Comparison between LSMM and GPA #####
# Vary K=10, 50, 100, 500,
# corr=0, 0.2, 0.4, 0.6 to get 
# Figures S56-S59 in Supplementary Document

library(LSMM)
library(pROC)
library(MASS)
library(GPA)

source("performance.R")

# function to generate data
generate_data_corr <- function(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, q){
  # design matrix of fixed effects
  Z         <- rep(0, M*L)
  indexZ    <- sample(M*L, M*L*Z.perc)
  Z[indexZ] <- 1
  Z         <- matrix(Z, M, L)
  
  # design matrix of random effects
  corr <- matrix(0, K, K)
  corr[1:10, 1:10] <- matrix(q, 10, 10)
  diag(corr) <- 1
  preA <- mvrnorm(M, rep(0, K), corr)
  A <- t(t(preA) < apply(preA, 2, quantile, probs = A.perc))
  
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
K      <- 500      # No. of random effects
Z.perc <- 0.1      # the proportion the entries in Z is 1
A.perc <- 0.1      # the proportion the entries in A is 1
alpha  <- 0.2      # parameter in the Beta distribution
beta0  <- -2       # intercept of the logistic model
set.seed(1)
b      <- rnorm(L) # fixed effects
omega  <- 0.2      # proportion of relevant annotations
sigma2 <- 1        # parameter in the spike-slab prior
corr   <- 0.2      # correlation coeffecient among the first 10 annotations with random effects
rep    <- 50       # repeat times

result <- matrix(0, rep, 12)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")
  
  data <- generate_data_corr(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, corr)
  
  # LSMM
  fit <- LSMM(data$Pvalue, data$Z, data$A)
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")
  result[i, 1:4]   <- as.numeric(performance(data$gamma, assoc.SNP$gamma, 1-fit$pi1))
  
  # LSMM without fixed effects
  fit1 <- LSMM(data$Pvalue, NULL, data$A)
  assoc1.SNP <- assoc.SNP(fit1, FDRset = 0.1, fdrControl = "global")
  result[i, 5:8]   <- as.numeric(performance(data$gamma, assoc1.SNP$gamma, 1-fit1$pi1))
  
  # GPA
  fit.GPA <- GPA(data$Pvalue, data$A)
  assoc.GPA <- assoc(fit.GPA, FDR = 0.1, fdrControl="global")
  result[i, 9:12] <- as.numeric(performance(data$gamma, assoc.GPA, fdr(fit.GPA)))
}

result <- as.data.frame(result)
names(result) <- c("FDR.LSMM", "power.LSMM", "AUC.LSMM", "pAUC.LSMM", 
                   "FDR.LSMM.no.fix", "power.no.fix", "AUC.no.fix", "pAUC.no.fix", 
                   "FDR.GPA", "power.GPA", "AUC.GPA", "pAUC.GPA")
