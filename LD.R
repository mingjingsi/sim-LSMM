##### Simulation study for evaluating the LD effects on LSMM #####
# Figure S44 in Supplementary Document

####### prepare ##########

# quality control
system("plink --bfile ctrl --geno 0.01 --mind 0.05 --hwe 0.001 --maf 0.05 --make-bed --out ctrlqc")

# use only SNPs in CHR 1
system("plink --bfile ctrlqc --chr 1 --make-bed --out ctrlqc_chr1")

# keep the SNPs with BP
SNP <- read.table(file = "ctrlqc_chr1.bim")
ref <- read.table(file = "wtccc.variant_function")
ref_chr1 <- ref[which(ref$V3 == "chr1"), ]
keep.SNP <- as.vector(SNP$V2[which(SNP$V2 %in% ref_chr1$V8)])
write.table(keep.SNP, file = "keep.SNP", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("plink --bfile ctrlqc_chr1 --extract keep.SNP --make-bed --out ctrl_qc_chr1")

SNP <- read.table(file = "ctrl_qc_chr1.bim")
write.csv(SNP, file = "ctrl_qc_chr1.csv")

write.csv(ref_chr1, file = "ref_chr1.csv")

# combine with Excel

SNP <- read.csv(file = "ctrl_qc_chr1.csv", header = FALSE)
write.table(SNP, file = "ctrl_qc_chr1.bim", quote = FALSE, row.names = FALSE, col.names = FALSE)

# set causal SNPs
SNP <- read.table(file = "ctrl_qc_chr1.bim")
set.seed(1)
causal.SNP <- sample(SNP$V2, round(length(SNP$V2)/1000))
write.table(causal.SNP, file = "causal.snplist", quote = FALSE, row.names = FALSE, col.names = FALSE)

# simulate GWAS data
system("gcta --bfile ctrl_qc_chr1  --simu-qt  --simu-causal-loci causal.snplist --simu-hsq 0.05 --simu-rep 50 --out sim_ctrl")

###### random annotation #####
library(LSMM)
tmp <- read.table(file = "tmp.fam")
sim <- read.table(file = "sim_ctrl.phen")
causal.SNP <- read.table(file = "causal.snplist")
SNP <- read.table(file = "tmp.bim")

causal.bp <- SNP$V4[which((SNP$V2 %in% causal.SNP$V1) == 1)]
distance <- 1e6
bp.range <- NULL
for (i in 1:length(causal.bp)){
  bp.range <- c(bp.range, seq(causal.bp[i] - distance, causal.bp[i] + distance))
}
index <- which((SNP$V4 %in% bp.range) == 1)

# sample Z & A
M <- nrow(SNP)
L <- 10
K <- 100
Z.perc <- 0.1
A.perc <- 0.1
Z         <- rep(0, M*L)
indexZ    <- sample(M*L, M*L*Z.perc)
Z[indexZ] <- 1
Z         <- matrix(Z, M, L)
Z.causal.perc <- 0.6
Z[index, 1:L] <- (matrix(runif(length(index)*L), length(index), L) < Z.causal.perc)
A         <- rep(0, M*K)
indexA    <- sample(M*K, M*K*A.perc)
A[indexA] <- 1
A         <- matrix(A, M, K)

est.random.K100 <- matrix(0, 50, M)
est.LFM.K100 <- matrix(0, 50, M)

for (i in 1:50){
  tmp$V6 <- sim[, i+2]
  write.table(tmp, file = "tmp.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp --assoc  --out sim_pvalue")
  pvalue <- read.table(file = "sim_pvalue.qassoc", header = TRUE)
  fit <- LSMM(pvalue$P, Z, A)
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")
  est.random.K100[i, ] <- assoc.SNP$gamma
  est.LFM.K100[i, ] <- assoc.SNP$gamma.stage2
}

save(est.random.K100, file = "est.random.K100.RData")
save(est.LFM.K100, file = "est.LFM.K100.RData")


##### annotation (within 1Mb) & no annotation #####
# SNPs around 1Mb
SNP <- read.table(file = "tmp.bim")
causal.SNP <- read.table(file = "causal.snplist")
causal.bp <- SNP$V4[which((SNP$V2 %in% causal.SNP$V1) == 1)]
distance <- 1e6
bp.range <- NULL
for (i in 1:length(causal.bp)){
  bp.range <- c(bp.range, seq(causal.bp[i] - distance, causal.bp[i] + distance))
}
index <- which((SNP$V4 %in% bp.range) == 1)

library(LSMM)
tmp <- read.table(file = "tmp.fam")
sim <- read.table(file = "sim_ctrl.phen")

M <- nrow(SNP)
L <- 10
K <- 100
Z.perc <- 0.1
A.perc <- 0.1
Z         <- rep(0, M*L)
indexZ    <- sample(M*L, M*L*Z.perc)
Z[indexZ] <- 1
Z         <- matrix(Z, M, L)
A         <- rep(0, M*K)
indexA    <- sample(M*K, M*K*A.perc)
A[indexA] <- 1
A         <- matrix(A, M, K)
Z.causal.perc <- 0.6
Z[index, 1:L] <- (matrix(runif(length(index)*L), length(index), L) < Z.causal.perc)
A.causal.perc <- 0.6
omega <- 0.2
A.col <- sample(K, K*omega)
A[index, A.col] <- (matrix(runif(length(index)*length(A.col)), length(index), length(A.col)) < A.causal.perc)

est.relevant.K100 <- matrix(0, 50, M)
est.woanno.K100 <- matrix(0, 50, M)

for (i in 1:50){
  tmp$V6 <- sim[, i+2]
  write.table(tmp, file = "tmp.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp --assoc  --out sim_pvalue")
  pvalue <- read.table(file = "sim_pvalue.qassoc", header = TRUE)
  fit <- LSMM(pvalue$P, Z, A)
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")
  est.relevant.K100[i, ] <- assoc.SNP$gamma
  est.woanno.K100[i, ] <- assoc.SNP$gamma.stage1
}

save(est.relevant.K100, file = "est.relevant.K100.RData")
save(est.woanno.K100, file = "est.woanno.K100.RData")

##### plot 15*15 #####
SNP <- read.table(file = "tmp.bim")
causal.SNP <- read.table(file = "causal.snplist")
causal.bp <- SNP$V4[which((SNP$V2 %in% causal.SNP$V1) == 1)]

FDR.summary <- NULL

methods <- c("woanno", "LFM", "random", "relevant")
Ks <- c("100", "500", "1000")
distance <- seq(1e5, 1e6, 1e5)

for(i in 1:length(methods)){
  for(j in 1:length(Ks)){
    load(paste("est.", methods[i], ".K", Ks[j], ".RData", sep = ""))
    est <- get(paste("est.", methods[i], ".K", Ks[j], sep = ""))
    FDR <- matrix(0, 50, 10)
    
    for (k in 1:length(distance)){
      bp.range <- NULL
      for (l in 1:length(causal.bp)){
        bp.range <- c(bp.range, seq(causal.bp[l] - distance[k], causal.bp[l] + distance[k]))
      }
      SNP.range <- SNP$V2[which((SNP$V4 %in% bp.range) == 1)]
      for (l in 1:50){
        TP <- sum(SNP$V2[which(est[l, ] == 1)] %in% SNP.range)
        
        if (sum(est[l, ]) == 0)
          FDR[l, k] <- 0
        else
          FDR[l, k] <- (sum(est[l, ]) - TP)/sum(est[l, ])
      }
    }
    
    FDR.summary$Method <- c(FDR.summary$Method, rep(methods[i], 50*10))
    FDR.summary$K <- c(FDR.summary$K, rep(Ks[j], 50*10))
    FDR.summary$distance <- c(FDR.summary$distance, rep(distance, each = 50))
    FDR.summary$FDR <- c(FDR.summary$FDR, as.vector(FDR))
    
  }
}

FDR.summary <- as.data.frame(FDR.summary)
FDR.summary$distance <- as.factor(FDR.summary$distance)
levels(FDR.summary$distance) <- c(paste(seq(100,1000,100), sep = ""))
levels(FDR.summary$K) <- c("K = 100", "K = 1000", "K = 500")
levels(FDR.summary$Method) <- c("Fixed effects", "Fixed + Random effects", 
                                "Fixed + Relevant random effects", "No effects")
FDR.summary$K <- factor(FDR.summary$K, levels = levels(FDR.summary$K)[c(1, 3, 2)])
FDR.summary$Method <- factor(FDR.summary$Method, levels = levels(FDR.summary$Method)[c(4, 1, 2, 3)])

save(FDR.summary, file = "FDR.summary.RData")

levels(FDR.summary$distance) <- c("100", " ", "300", "  ", "500", "   ", "700", "    ", "900",
                                  "     ")

library(ggplot2)

p <- ggplot(FDR.summary, aes(x = distance, y = FDR)) +
  geom_boxplot() +
  facet_grid(Method ~ K) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
  ylim(0, 0.5) +
  xlab("Distance (kb)") +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 19)) 
