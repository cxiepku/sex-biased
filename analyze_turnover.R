library(matrixStats)    # 1.3.0
library(dplyr)    # 1.1.4


# function to calculate regulatory turnover of all genes between taxa
calculate.turnover <- function(TPM1.df, TPM2.df, coldata1.df, coldata2.df, gene.metadata.cols,
                               filter.cutoff, biased.cutoff, BH.cutoff, out.name) {
    GeneID.df <- data.frame(GeneID = rownames(TPM1.df))
    result.df <- cbind(GeneID.df, TPM1.df[, gene.metadata.cols])

    for (sex in c('F', 'M')) {
        coldata.df <- rbind(coldata1.df, coldata2.df) %>% filter(Sex == sex)

## add 1 on TPM value
        TPM.plus1 <- as.matrix(cbind(TPM1.df, TPM2.df)[, rownames(coldata.df)] + 1)
        
        tmp.df <- data.frame(GeneID = rownames(TPM1.df))

## calculate medians and ratio
        tmp.df$Median_s1 <- rowMedians(TPM.plus1[, rownames(coldata.df[coldata.df$Species == species1, ])])
        tmp.df$Median_s2 <- rowMedians(TPM.plus1[, rownames(coldata.df[coldata.df$Species == species2, ])])
        tmp.df$Ratio <- tmp.df$Median_s1 / tmp.df$Median_s2

## perform Wilcoxon rank sum test
        tmp.df$P <- sapply(1:nrow(TPM.plus1), function(i) {
            row.data <- cbind.data.frame(Gene = as.numeric(t(TPM.plus1[i, ])), Condition = coldata.df$Species)
            row.p <- wilcox.test(Gene ~ Condition, row.data)$p.value
            return(row.p)
        })
        
## identify turnover genes
        tmp.df$P_adj <- NA
        
        filtered.index <- (tmp.df$Median_s1 > filter.cutoff) | (tmp.df$Median_s2 > filter.cutoff)
        turnover.index <- filtered.index & ((tmp.df$Ratio > biased.cutoff) | (tmp.df$Ratio < 1 / biased.cutoff))
        
        tmp.df[turnover.index, 'P_adj'] <- p.adjust(tmp.df[turnover.index, 'P'], method = 'BH')    ## FDR
        turnover.index <- turnover.index & !(is.na(tmp.df$P_adj)) & (tmp.df$P_adj < BH.cutoff)
        
        if (sex == 'F') {
            result.df$F_filtered <- filtered.index
            result.df$F_turnover <- turnover.index
        } else {
            result.df$M_filtered <- filtered.index
            result.df$M_turnover <- turnover.index
        }
    }
    
    write.table(result.df, out.name, quote = F, sep = '\t', row.names = F)
}


# input parameters
gene.metadata.cols <- 1:7    ## gene metadata columns

filter.cutoff <- 2    ## median TPM value in at least one taxon is larger than 1
biased.cutoff <- 1.25    ## ratio cutoff for turnover genes
BH.cutoff <- 0.1    ## FDR (BH) cutoff


species1 <- 'DOM'
species2 <- 'MUS'
out.name <- paste0('data/', species1, '_', species2, '_brain_turnover.tsv')    ## output file name

# input files
## TPM files with gene metadata
TPM1.df <- read.table(paste0('data/', species1, '_brain_TPM.tsv'), header = T, sep = '\t', row.names = 'GeneID')
TPM2.df <- read.table(paste0('data/', species2, '_brain_TPM.tsv'), header = T, sep = '\t', row.names = 'GeneID')

## sample annotation files
coldata1.df <- read.table(paste0('data/', species1, '_brain_coldata.tsv'), header = T, sep = '\t', row.names = 'FN')
coldata2.df <- read.table(paste0('data/', species2, '_brain_coldata.tsv'), header = T, sep = '\t', row.names = 'FN')

## identify regulatory turnover of all genes between taxa
calculate.turnover(TPM1.df, TPM2.df, coldata1.df, coldata2.df, gene.metadata.cols,
                   filter.cutoff, biased.cutoff, BH.cutoff, out.name)

