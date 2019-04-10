##### Performance of LSMM when the proportion of risk SNPs was extremely small #####
# Vary alpha=0.2, 0.4, 0.6, 
# pi1=0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2 to get 
# Figures S45 in Supplementary Document

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

# Higher Criticism
epsest.func <- function(x,u,sigma)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation
  
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)
  
  epsest=NULL
  
  for (j in 1:length(tt)) { 
    
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi
    
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    } 
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}

M     <- 100000   # No. of SNPs
alpha <- 0.2      # parameter in the Beta distribution
pi1   <- 0.1      # proportion of risk SNPs
rep   <- 50       # repeat times

result <- matrix(0, rep, 3)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")
  
  data <- generate_data_TGM(M, alpha, pi1)

  # default
  fit1 <- LSMM(data$Pvalue)
  assoc.SNP1 <- assoc.SNP(fit1, FDRset = 0.1, fdrControl = "global")
  
  result[i, 1] <- fit1$pi1_
  result[i, 2] <- perf_FDR(data$gamma, assoc.SNP1$gamma, 1-fit1$pi1)
  
  # HC
  zscore <- (qnorm(data$Pvalue/2))*(-1)^(runif(M)<0.5)

  result[i, 3] <- epsest.func(zscore, 0, 1)
  
}

result <- as.data.frame(result)
names(result) <- c("est.LSMM", "FDR.LSMM", "est.HC")
