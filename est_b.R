##### Estimation of b #####
# Vary alpha=0.2, 0.4, 0.6, 
# K=100, 500, 1000,
# omega=0, 0.25, 0.5, 0.75, 1 to get 
# Figures S28-S38 in Supplementary Document

library(LSMM)
library(MASS)

# function to generate data
generate_data <- function(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2){
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
omega  <- 0.25     # proportion of relevant annotations
sigma2 <- 1        # parameter in the spike-slab prior
rep    <- 50       # repeat times

result1 <- matrix(0, rep, 11)
result2 <- matrix(0, rep, 11)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")
  
  data <- generate_data(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2)
  
  fit.LSMM <- LSMM(data$Pvalue, data$Z, data$A)
  fit.LFM  <- LSMM(data$Pvalue, data$Z, NULL)
  
  result1[i, 1:11] <- as.numeric(fit.LSMM$b)
  result2[i, 1:11] <- as.numeric(fit.LFM$b)
}

result1 <- as.data.frame(result1)
names(result1) <- c("b0.LSMM", "b1.LSMM", "b2.LSMM", "b3.LSMM", "b4.LSMM", "b5.LSMM",
                    "b6.LSMM", "b7.LSMM", "b8.LSMM", "b9.LSMM", "b10.LSMM")

result2 <- as.data.frame(result2)
names(result2) <- c("b0.LFM", "b1.LFM", "b2.LFM", "b3.LFM", "b4.LFM", "b5.LFM",
                    "b6.LFM", "b7.LFM", "b8.LFM", "b9.LFM", "b10.LFM")

