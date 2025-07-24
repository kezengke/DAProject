# Plot log10 pvalue plots
rm(list = ls())
library(MetagenomeTools)
library(ggplot2)
library(dplyr)
library(patchwork)

# RDP
all_names<-gsub(basename(list.files("CountsTables/RDPRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/RDP/Log10PvPPlots(", name, ").png"), width=5000, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/RDP/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/RDP/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/RDP/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/RDP/edgeR/", name, "_edger.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)

  # tvw
  p1<-Log10PvPPlot(tRES, wRES, "t-test", "Wilcoxon", paste0("(RDP-", name, ") t-test vs. Wilcoxon"))
  # tvd
  p2<-Log10PvPPlot(tRES, dRES, "t-test", "DESeq2", paste0("(RDP-", name, ") t-test vs. DESeq2"))
  # tve
  p3<-Log10PvPPlot(tRES, eRES, "t-test", "edgeR", paste0("(RDP-", name, ") t-test vs. edgeR"))

  # wvd
  p4<-Log10PvPPlot(wRES, dRES, "Wilcoxon", "DESeq2", paste0("(RDP-", name, ") Wilcoxon vs. DESeq2"))
  # wve
  p5<-Log10PvPPlot(wRES, eRES, "Wilcoxon", "edgeR", paste0("(RDP-", name, ") Wilcoxon vs. edgeR"))

  # dve
  p6<-Log10PvPPlot(dRES, eRES, "DESeq2", "edgeR", paste0("(RDP-", name, ") DESeq2 vs. edgeR"))

  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}


# DADA2
all_names<-gsub(basename(list.files("CountsTables/dada2Raw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/dada2/Log10PvPPlots(", name, ").png"), width=5000, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/dada2/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/dada2/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/dada2/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/dada2/edgeR/", name, "_edger.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)

  # tvw
  p1<-Log10PvPPlot(tRES, wRES, "t-test", "Wilcoxon", paste0("(dada2-", name, ") t-test vs. Wilcoxon"))
  # tvd
  p2<-Log10PvPPlot(tRES, dRES, "t-test", "DESeq2", paste0("(dada2-", name, ") t-test vs. DESeq2"))
  # tve
  p3<-Log10PvPPlot(tRES, eRES, "t-test", "edgeR", paste0("(dada2-", name, ") t-test vs. edgeR"))

  # wvd
  p4<-Log10PvPPlot(wRES, dRES, "Wilcoxon", "DESeq2", paste0("(dada2-", name, ") Wilcoxon vs. DESeq2"))
  # wve
  p5<-Log10PvPPlot(wRES, eRES, "Wilcoxon", "edgeR", paste0("(dada2-", name, ") Wilcoxon vs. edgeR"))

  # dve
  p6<-Log10PvPPlot(dRES, eRES, "DESeq2", "edgeR", paste0("(dada2-", name, ") DESeq2 vs. edgeR"))

  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}

# WGS
all_names<-gsub(basename(list.files("CountsTables/WGSRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/WGS/Log10PvPPlots(", name, ").png"), width=5000, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/WGS/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/WGS/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/WGS/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/WGS/edgeR/", name, "_edger.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)

  # tvw
  p1<-Log10PvPPlot(tRES, wRES, "t-test", "Wilcoxon", paste0("(WGS-", name, ") t-test vs. Wilcoxon"))
  # tvd
  p2<-Log10PvPPlot(tRES, dRES, "t-test", "DESeq2", paste0("(WGS-", name, ") t-test vs. DESeq2"))
  # tve
  p3<-Log10PvPPlot(tRES, eRES, "t-test", "edgeR", paste0("(WGS-", name, ") t-test vs. edgeR"))

  # wvd
  p4<-Log10PvPPlot(wRES, dRES, "Wilcoxon", "DESeq2", paste0("(WGS-", name, ") Wilcoxon vs. DESeq2"))
  # wve
  p5<-Log10PvPPlot(wRES, eRES, "Wilcoxon", "edgeR", paste0("(WGS-", name, ") Wilcoxon vs. edgeR"))

  # dve
  p6<-Log10PvPPlot(dRES, eRES, "DESeq2", "edgeR", paste0("(WGS-", name, ") DESeq2 vs. edgeR"))

  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}

# RNAseq
all_names<-gsub(basename(list.files("CountsTables/RNAseqRaw/", pattern = "*.txt",full.names = TRUE)), pattern=".txt$", replacement="")

for(name in all_names){
  png(paste0("Plots/RNAseq/Log10PvPPlots(", name, ").png"), width=5000, height=2400, res = 300)
  par(mar=c(5,6,4,1)+.1)
  tRES<-read.table(paste0("PkgResults/RNAseq/ttest/", name, "_t.txt"), header = T, row.names = 1)
  wRES<-read.table(paste0("PkgResults/RNAseq/Wilcoxon/", name, "_wilcox.txt"), header = T, row.names = 1)
  dRES<-read.table(paste0("PkgResults/RNAseq/DESeq2/", name, "_deseq2.txt"), header = T, row.names = 1)
  eRES<-read.table(paste0("PkgResults/RNAseq/edgeR/", name, "_edger.txt"), header = T, row.names = 1)

  tRES<-processTtestRes(tRES)
  wRES<-processWilcoxonRes(wRES)
  dRES<-processDESeq2Res(dRES)
  eRES<-processEdgeRRes(eRES)

  # tvw
  p1<-Log10PvPPlot(tRES, wRES, "t-test", "Wilcoxon", paste0("(RNAseq-", name, ") t-test vs. Wilcoxon"))
  # tvd
  p2<-Log10PvPPlot(tRES, dRES, "t-test", "DESeq2", paste0("(RNAseq-", name, ") t-test vs. DESeq2"))
  # tve
  p3<-Log10PvPPlot(tRES, eRES, "t-test", "edgeR", paste0("(RNAseq-", name, ") t-test vs. edgeR"))

  # wvd
  p4<-Log10PvPPlot(wRES, dRES, "Wilcoxon", "DESeq2", paste0("(RNAseq-", name, ") Wilcoxon vs. DESeq2"))
  # wve
  p5<-Log10PvPPlot(wRES, eRES, "Wilcoxon", "edgeR", paste0("(RNAseq-", name, ") Wilcoxon vs. edgeR"))

  # dve
  p6<-Log10PvPPlot(dRES, eRES, "DESeq2", "edgeR", paste0("(RNAseq-", name, ") DESeq2 vs. edgeR"))

  combined_plots <- p1 + p2 + p3 + p4 + p5 + p6
  print(combined_plots + plot_layout(ncol = 3))

  dev.off()
}
