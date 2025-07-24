# Calculating wilcoxon results
rm(list = ls())
source("CalculationTools.R")
library(coin)

# RDP
rm(list = ls())
all_files<-list.files("CountsTables/RDPNorm", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  results<-calcWilcox(countsT, meta)

  write.table(results, paste0("PkgResults/RDP/Wilcoxon/", gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
}

# dada2
rm(list = ls())
all_files<-list.files("CountsTables/dada2Norm", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  results<-calcWilcox(countsT, meta)

  write.table(results, paste0("PkgResults/dada2/Wilcoxon/", gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
}

# WGS
rm(list = ls())
all_files<-list.files("CountsTables/WGSNorm", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  results<-calcWilcox(countsT, meta)

  write.table(results, paste0("PkgResults/WGS/Wilcoxon/", gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
}

# RNAseq
rm(list = ls())
all_files<-list.files("CountsTables/RNAseqNorm", pattern = "*.txt",full.names = TRUE)

for (file in all_files) {
  countsT<-read.table(file, sep = "\t", header = T, row.names = 1, check.names = F)
  meta<-read.table(paste0("MetaData/metadata_", gsub(basename(file), pattern=".txt$", replacement=""), ".txt"),
                   header = T, row.names = 1)

  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]

  rownames(meta)<-colnames(countsT)
  colnames(meta)<-"conditions"

  results<-calcWilcox(countsT, meta)

  write.table(results, paste0("PkgResults/RNAseq/Wilcoxon/", gsub(basename(file), pattern=".txt$", replacement=""), "_wilcox.txt"), sep = "\t", row.names = T)
}


