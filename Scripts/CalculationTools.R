# Load required packages quietly
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(coin)
})

#' Load in counts table
LoadCountsT <- function(filePath) {
  if (grepl("\\.txt$", filePath, ignore.case = TRUE)) {
    countsT<-read.table(filePath, sep = "\t", header = T, row.names = 1, check.names = F)
    return(countsT)
  } else if (grepl("\\.csv$", filename, ignore.case = TRUE)) {
    countsT<-read.csv(file, header = T, row.names = 1, check.names = F)
    return(countsT)
  } else {
    stop("unsupported filetype")
  }
}

#' Load in metadata
LoadMeta <- function(filePath) {
  if (grepl("\\.txt$", filePath, ignore.case = TRUE)) {
    meta<-read.table(filePath, sep = "\t", header = T, row.names = 1, check.names = F)
    return(meta)
  } else if (grepl("\\.csv$", filename, ignore.case = TRUE)) {
    meta<-read.csv(file, header = T, row.names = 1, check.names = F)
    return(meta)
  } else {
    stop("unsupported filetype")
  }
}

#' Function to normalize counts table
normFun <- function(table) {
  n<-colSums(table)
  sumx<-sum(table)
  for (j in 1:ncol(table)) {
    table[,j]<-table[,j]/n[j]
  }
  table<-log10(table*(sumx/ncol(table))+1)
  table<-data.frame(table, check.names = F)
  return(table)
}

calcTtest <- function(table, meta) {
  t_stats<-apply(table, 1, function(x){t.test(unlist(x)~meta$conditions)$stat})
  t_test_p<-apply(table, 1, function(x){t.test(unlist(x)~meta$conditions)$p.value})

  t_results<-cbind(t_stats, t_test_p)
  rownames(t_results)<-rownames(table)
  colnames(t_results)<-c("stats", "pval")
  t_results<-data.frame(t_results, check.names = F)
  return(t_results)
}

#' Function to calculate Wilcoxon results
calcWilcox <- function(table, meta) {
  wilcox_stats<-apply(table, 1, function(x){statistic(wilcox_test(unlist(x)~factor(meta$conditions)))})
  wilcox_p<-apply(table, 1, function(x){pvalue(wilcox_test(unlist(x)~factor(meta$conditions)))})

  wilcox_results<-cbind(wilcox_stats, wilcox_p)
  rownames(wilcox_results)<-rownames(table)
  colnames(wilcox_results)<-c("stats", "pval")
  wilcox_results<-data.frame(wilcox_results, check.names = F)
  return(wilcox_results)
}

#' Function to calculate DESeq2 results
calcDESeq2 <- function(table, meta) {
  #solve deseq2 all 0 issue
  table<-table+1

  meta$conditions<-factor(meta$conditions)
  dds1 <- DESeqDataSetFromMatrix(countData=table,
                                 colData=meta,
                                 design=~conditions)

  # Wrap DESeq in tryCatch to handle errors
  dds2 <- tryCatch({
    DESeq(dds1)
  }, error = function(e) {
    # If error occurs, return NA instead of running DESeq
    return(NULL)
  })

  # If dds2 is NULL, skip the remaining lines and create deseq_results as NA
  if (is.null(dds2)) {
    deseq_results <- matrix(NA, nrow = nrow(table), ncol = 2)
  } else {
    res <- results(dds2, cooksCutoff=FALSE, independentFiltering=FALSE)
    deseq_results<-cbind(res$stat, res$pvalue)
  }

  rownames(deseq_results)<-rownames(table)
  colnames(deseq_results)<-c("stats", "pval")
  deseq_results<-data.frame(deseq_results, check.names = F)
  return(deseq_results)
}

#' Function to calculate edgeR results
calcEdgeR <- function(table, meta) {
  group <- meta$condition
  dgList <- DGEList(counts=table, group = group)
  dgList <- calcNormFactors(dgList, method="TMM")

  # Wrap estimateDisp in tryCatch to handle errors
  dgList <- tryCatch({
    estimateDisp(dgList)
  }, error = function(e) {
    # If error occurs, return NA instead of running estimateDisp
    return(NULL)
  })

  # If dgList is NULL, skip the remaining lines and create edger_results as NA
  if (is.null(dgList)) {
    edger_results <- matrix(NA, nrow = nrow(table), ncol = 2)
  } else {
    et <- exactTest(dgList)
    res <- et$table
    edger_results <- cbind(res$logFC, res$PValue)
  }

  rownames(edger_results) <- rownames(table)
  colnames(edger_results) <- c("stats", "pval")
  edger_results <- data.frame(edger_results, check.names = F)

  return(edger_results)
}

#' Function to resample counts table (whole taxon) with multiple (Rnbinom)
resampleWholeTaxonRNBINOM <- function(table, meta, multiple) {
  if (nrow(table) == 0 || nrow(meta) == 0) {
    stop("Input table or meta data frame is empty.")
  }

  # Function to calculate mean and sd for each group
  calculateMeanVar <- function(z) {
    c(meanTotal = mean(z), varTotal = var(z))
  }

  # Apply the function to each row of the table and combine results into a data frame
  MeanVar_table <- t(apply(table, 1, calculateMeanVar))
  MeanVar_table <- as.data.frame(MeanVar_table)
  rownames(MeanVar_table) <- rownames(table)

  exclude<-which(MeanVar_table$varTotal <= MeanVar_table$meanTotal)
  # Function to generate resampled counts for each row
  resample_counts <- function(row_index) {
    z <- table[row_index,]

    n<-length(z)

    r<-((MeanVar_table[row_index, 1])^2)/((MeanVar_table[row_index, 2])-MeanVar_table[row_index, 1])
    p<-r/(r+(MeanVar_table[row_index, 1]))

    nz<-rnbinom(n = n, size = r,  prob = p)

    return(nz)
  }

  if(length(exclude)>0){
    table<-table[-exclude, , drop = F]
    MeanVar_table<-MeanVar_table[-exclude, , drop = F]
  }

  # Apply the resampling function to each row
  newT <- t(sapply(seq_len(nrow(table)), resample_counts))

  # Set column and row names
  colnames(newT) <- colnames(table)
  newT <- newT[, colnames(table)]
  rownames(newT) <- rownames(table)

  # # Replace negative values with 0 and round to integer
  # newT[newT < 0] <- 0
  # Add smallest number to whole counts table
  # newT<-newT + abs(min(newT))
  newT <- data.frame(round(newT, digits = 0), check.names = F)

  return(newT)
}
