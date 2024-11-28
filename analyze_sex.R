library(matrixStats)    # 1.3.0
library(dplyr)    # 1.1.4


# function to identify and analyze sex-biased genes
analyze.sex <- function(TPM.df, coldata.df, gene.metadata.cols, excluded.chroms,
                        filter.cutoff, biased.cutoff, unbiased.cutoff, BH.cutoff, out.prefix) {
## add 1 on TPM value
    TPM.plus1 <- as.matrix(TPM.df[, rownames(coldata.df)] + 1)

    TPM.plus1.F <- TPM.plus1[, rownames(coldata.df[coldata.df$Sex == 'F', ])]
    TPM.plus1.M <- TPM.plus1[, rownames(coldata.df[coldata.df$Sex == 'M', ])]
    
    GeneID.df <- data.frame(GeneID = rownames(TPM.df))
    result.df <- cbind(GeneID.df, TPM.df[, gene.metadata.cols])

## calculate Median_F, Median_M, and Ratio
    result.df$Median_F <- rowMedians(TPM.plus1.F)
    result.df$Median_M <- rowMedians(TPM.plus1.M)
    result.df$Ratio <- result.df$Median_F / result.df$Median_M
    
## perform Wilcoxon rank sum test
    result.df$P <- sapply(1:nrow(TPM.plus1), function(i) {
        row.data <- cbind.data.frame(Gene = as.numeric(t(TPM.plus1[i, ])), Condition = coldata.df$Sex)
        row.p <- wilcox.test(Gene ~ Condition, row.data)$p.value
        return(row.p)
    })

## identify sex-biased and unbiased genes
    result.df$P_adj <- NA
    
    filtered.index <- (result.df$Median_F > filter.cutoff) | (result.df$Median_M > filter.cutoff)
    F.index <- filtered.index & (result.df$Ratio > biased.cutoff)
    M.index <- filtered.index & (result.df$Ratio < 1 / biased.cutoff)
    U.index <- filtered.index & (result.df$Ratio < unbiased.cutoff) & (result.df$Ratio > 1 / unbiased.cutoff)
    
    result.df[F.index | M.index, 'P_adj'] <- p.adjust(result.df[F.index | M.index, 'P'], method = 'BH')    ## FDR
    F.index <- F.index & !(is.na(result.df$P_adj)) & (result.df$P_adj < BH.cutoff)
    M.index <- M.index & !(is.na(result.df$P_adj)) & (result.df$P_adj < BH.cutoff)

    result.df$Type <- NA
    result.df[F.index, 'Type'] <- 'F'    ## female-biased
    result.df[M.index, 'Type'] <- 'M'    ## male-biased
    result.df[U.index, 'Type'] <- 'U'    ## unbiased

## calculate IQR_F, IQR_M, and Disp = IQR / Median
    result.df$IQR_F <- rowIQRs(TPM.plus1.F)
    result.df$IQR_M <- rowIQRs(TPM.plus1.M)

    F.disp <- result.df$IQR_F / result.df$Median_F
    M.disp <- result.df$IQR_M / result.df$Median_M

    result.df$Disp <- NA
    result.df[F.index, 'Disp'] <- F.disp[F.index]
    result.df[M.index, 'Disp'] <- M.disp[M.index]
    result.df[U.index, 'Disp'] <- (F.disp[U.index] + M.disp[U.index]) / 2

    write.table(result.df, paste0(out.prefix, 'sex.tsv'), quote = F, sep = '\t', row.names = F)

## calculate sex-bias index (SBI)
    M.index.woY <- M.index & !(result.df$Chr %in% excluded.chroms)

    if (sum(F.index) > 0 & sum(M.index.woY) > 0) {
        coldata.df$FN <- rownames(coldata.df)
        coldata.df$SBI <- colMedians(TPM.plus1[F.index, ]) - colMedians(TPM.plus1[M.index.woY, ])
        coldata.df$SBI_scaled <- 2 * (coldata.df$SBI - min(coldata.df$SBI)) / (max(coldata.df$SBI) - min(coldata.df$SBI)) - 1
        
        write.table(coldata.df, paste0(out.prefix, 'SBI.tsv'), quote = F, sep = '\t', row.names = F)
    }
}


# input parameters
gene.metadata.cols <- 1:7    ## gene metadata columns
excluded.chroms <- c('Y')    ## chromosomes excluded for sex-bias index (SBI) calculation

filter.cutoff <- 2    ## median TPM value in at least one sex is larger than 1
biased.cutoff <- 1.25    ## ratio cutoff for sex-biased genes
unbiased.cutoff <- 1.05    ## ratio cutoff for unbiased genes
BH.cutoff <- 0.1    ## FDR (BH) cutoff


for (species in c('DOM', 'MUS')) {
    out.prefix <- paste0('data/', species, '_brain_')    ## prefix of output file names
    
# input files
## TPM file with gene metadata
    TPM.df <- read.table(paste0('data/', species, '_brain_TPM.tsv'), header = T, sep = '\t', row.names = 'GeneID')

## sample annotation file
    coldata.df <- read.table(paste0('data/', species, '_brain_coldata.tsv'), header = T, sep = '\t', row.names = 'FN')
    coldata.df$Sex <- factor(coldata.df$Sex, levels = c('M', 'F'))

## identify and analyze sex-biased genes
    analyze.sex(TPM.df, coldata.df, gene.metadata.cols, excluded.chroms,
                filter.cutoff, biased.cutoff, unbiased.cutoff, BH.cutoff, out.prefix)
}

