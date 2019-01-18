##### An example to reproduce the results of real data analysis #####
# Schizophrenia 2 (Schizophrenia Psychiatric Genome-Wide Association Study (GWAS) Consortium, 2011)
# For other GWASs, the code is similar only with some modifications explained below
# The data sets analyzed in this paper, including nine genic category annotations, 127 cell-type specific functional annotations, and the summary statistics of 30 GWAS can be downloaded at https://drive.google.com/drive/folders/1YrhncaNpo6doHfnIYqXPlhPdE7YpqEjO

library(LSMM)

load("SCZ1.RData") # load the summary statistics (notes: change the file name for other GWASs)
load("region9.RData") # load the genic category annotations
load("tissue127.RData") # load the cell-type specific annotations

SCZ1_Pvalue <- SCZ1[, c(1, 8)] # extract the SNP IDs and P-values (notes: the index of the columns may be different for other GWASs)

region9tissue127 <- as.data.frame(cbind(region9, tissue127))
region9tissue127$SNP <- rownames(region9)
  
SCZ1_ZA1 <- merge(SCZ1_Pvalue, region9tissue127, by.x = "snpid", by.y = "SNP") # keep the SNPs shared in GWAS and annotations (notes: the column name of the SNP ID may be different for other GWASs)

Pvalue <- SCZ1_ZA1[, 2]
Z <- SCZ1_ZA1[, 3:11]
A <- SCZ1_ZA1[, 12:138]

fit.LSMM <- LSMM(Pvalue, Z, A)

assoc.SNP.LSMM <- assoc.SNP(fit.LSMM, FDRset = 0.1, fdrControl="global")
table(assoc.SNP.LSMM$gamma)

relev.Anno.LSMM <- relev.Anno(fit.LSMM, FDRset = 0.1, fdrControl="local")
table(relev.Anno.LSMM)
