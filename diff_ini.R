##### Simulations to investigate the sensitivity of LSMM to initial parameter specification #####
# Vary alpha=0.2, 0.4, 0.6, 
# pi1=0.01, 0.05, 0.1, 0.15, 0.2 to get 
# Figures S24 in Supplementary Document

library(LSMM)
library(MASS)

# function to generate data
generate_data_TGM <- function(M, alpha, pi1){
  index <- sample(M, M*pi1)
  
  gamma <- numeric(M)
  gamma[index] <- 1
  
  Pvalue        <- runif(M)
  Pvalue[index] <- rbeta(length(index), alpha, 1)
  
  return(list(Pvalue = Pvalue, gamma = gamma))
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

M     <- 100000   # No. of SNPs
alpha <- 0.2      # parameter in the Beta distribution
pi1   <- 0.1      # proportion of risk SNPs
rep   <- 50       # repeat times

result <- matrix(0, rep, 4)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")
  
  data <- generate_data_TGM(M, alpha, pi1)

  # default
  fit1 <- LSMM(data$Pvalue)
  assoc.SNP1 <- assoc.SNP(fit1, FDRset = 0.1, fdrControl = "global")
  
  result[i, 1] <- fit1$pi1_
  result[i, 2] <- perf_FDR(data$gamma, assoc.SNP1$gamma, 1-fit1$pi1)
  
  # random initial values
  alpha_ini <- runif(1, 0.1, 0.6)
  pi1_ini <- runif(1, 0, 0.3)
  fit2 <- LSMM(data$Pvalue, alpha = alpha_ini, pi1_ = pi1_ini)
  assoc.SNP2 <- assoc.SNP(fit2, FDRset = 0.1, fdrControl = "global")
  
  result[i, 3] <- fit2$pi1_
  result[i, 4] <- perf_FDR(data$gamma, assoc.SNP2$gamma, 1-fit2$pi1)
  
}

result <- as.data.frame(result)
names(result) <- c("est.default", "FDR.default", "est.random", "FDR.random")
