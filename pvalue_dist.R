##### Simulations if p-values are not from beta distribution #####
# Vary K=100, 500, 1000,
# r=4, 1, 1/4,
# dist="spiky", "near_normal", "skew", "big_normal" to get
# Figures S41-S43 in Supplementary Document

library(LSMM)
library(MASS)

# function to generate data
generate_data_probit_dist <- function(M, L, K, Z.perc, A.perc, beta0, b, omega, sigma2, r, dist){
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
  sigmae2           <- var(Z %*% b + A %*% beta)/r # r is signal-noise ratio
  y                 <- beta0 + Z %*% b + A %*% beta + sqrt(sigmae2) * rnorm(M)
  gamma             <- rep(0, M)
  indexgamma        <- (y > 0)
  gamma[indexgamma] <- 1
  
  # Pvalue
  Pvalue <- runif(M)
  fz <- match.fun(dist)
  z <- fz(sum(indexgamma))
  Pvalue[indexgamma] <- pnorm(abs(z), lower.tail = FALSE)*2
  
  return( list(Z = Z, A = A, Pvalue = Pvalue, beta = beta, eta = eta, gamma = gamma))
}

spiky <- function(N){
  r <- runif(N)
  z <- rnorm(N, 0, 0.25)
  if(sum(r <= 0.2) != 0)
    z[which(r <= 0.2)] <- rnorm(sum(r <= 0.2), 0, 0.5)
  if(sum(r > 0.2 & r <= 0.4) != 0)
    z[which(r > 0.2 & r <= 0.4)] <- rnorm(sum(r > 0.2 & r <= 0.4), 0, 1)
  if(sum(r > 0.8) != 0)
    z[which(r > 0.8)] <- rnorm(sum(r > 0.8), 0, 2)
  
  return(z)
}

near_normal <- function(N){
  r <- runif(N)
  z <- rnorm(N, 0, 1)
  if(sum(r <= 1/3) != 0)
    z[which(r <= 1/3)] <- rnorm(sum(r <= 1/3), 0, 2)
  
  return(z)
}

skew <- function(N){
  r <- runif(N)
  z <- rnorm(N, -2, 2)
  if(sum(r <= 0.25) != 0)
    z[which(r <= 0.25)] <- rnorm(sum(r <= 0.25), -1, 1.5)
  if(sum(r > 0.25 & r <= (0.25+1/3)) != 0)
    z[which(r > 0.25 & r <= (0.25+1/3))] <- rnorm(sum(r > 0.25 & r <= (0.25+1/3)), 0, 1)
  if(sum(r > 5/6) != 0)
    z[which(r > 5/6)] <- rnorm(sum(r > 5/6), 1, 1)
  
  return(z)
}

big_normal <- function(N){
  r <- runif(N)
  z <- rnorm(N, 0, 4)
  
  return(z)
}

# evaluate FDR
perf_FDR <- function(true, est, fdr){
  
  fdr <- as.numeric(fdr)
  t <- table(true, est)
  if (sum(est)==0){
    FDR.fit <- 0
  }
  else if (sum(est)==length(est)){
    FDR.fit <- t[1]/(t[1]+t[2])
  }
  else{
    FDR.fit <- t[1,2]/(t[1,2]+t[2,2])
  }
  
  return(FDR.fit)
}

M      <- 100000   # No. of SNPs
L      <- 10       # No. of fixed effects
K      <- 100      # No. of random effects
Z.perc <- 0.1      # the proportion the entries in Z is 1
A.perc <- 0.1      # the proportion the entries in A is 1
alpha  <- 0.2      # parameter in the Beta distribution
beta0  <- -1       # intercept of the probit model
set.seed(1)
b      <- rnorm(L) # fixed effects
omega  <- 0.2      # proportion of relevant annotations
sigma2 <- 1        # parameter in the spike-slab prior
r      <- 4        # signal-noise level of probit model
dist   <- "spiky"  # distributions for z-scores
rep    <- 50       # repeat times

result <- matrix(0, rep, 3)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")
  
  data <- generate_data_probit_dist(M, L, K, Z.perc, A.perc, beta0, b, omega, sigma2, r, dist)
  
  fit <- LSMM(data$Pvalue, data$Z, data$A)
  
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")
  result[i, 1] <- perf_FDR(data$gamma, assoc.SNP1$gamma, 1-fit1$pi1)
  result[i, 2] <- perf_FDR(data$gamma, assoc.SNP$gamma.stage1, 1-fit$pi1.stage1)
  result[i, 3] <- perf_FDR(data$gamma, assoc.SNP$gamma.stage2, 1-fit$pi1.stage2)
  
}

result <- as.data.frame(result)
names(result) <- c("FDR.LSMM", "FDR.TGM", "FDR.LFM")
