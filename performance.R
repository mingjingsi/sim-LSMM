##### Function to measure performance #####
library(pROC)

performance <- function(true, est, fdr){
  
  fdr <- as.numeric(fdr)
  t <- table(true, est)
  if (sum(est)==0){
    FDR.fit <- 0
    power   <- 0
  }
  else if (sum(est)==length(est)){
    FDR.fit <- t[1]/(t[1]+t[2])
    power   <- 1
  }
  else{
    FDR.fit <- t[1,2]/(t[1,2]+t[2,2])
    power   <- t[2,2]/(t[2,1]+t[2,2])
  }
  AUC <- as.numeric(roc(true, fdr)$auc)
  pAUC <- as.numeric(roc(true, fdr, partial.auc = c(1, 0.8), parial.auc.correct = TRUE)$auc)
  
  return( list( FDR.fit = FDR.fit, power = power, AUC = AUC, pAUC = pAUC))
}

