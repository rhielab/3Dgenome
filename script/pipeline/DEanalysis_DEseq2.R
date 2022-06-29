library(DESeq2)
library(dplyr)
final <- read.csv(file = '/scratch/beoungle/probing_rnaseq/c42b/output/c42b/HTseq-GenesCountsTable.csv', header = T,check.names = F, row.names = 1)
expdesign <- read.csv(file = '/scratch/beoungle/probing_rnaseq/c42b/output/c42b/c42b_MIX_expdesign.csv', header = T, check.names = F, row.names = 1)
expdesign <- expdesign[match(names(final), rownames(expdesign)),] %>% as.data.frame() 
rownames(expdesign) <- names(final) 

names(expdesign) <- 'condition' 
condition <- levels(factor(expdesign$condition)) 
condition.com <- combn(condition,2) 
if(any(grepl('control',condition.com,ignore.case = T)) == T) {
condition.com <- condition.com[,apply(condition.com, 2, function(x) any(grepl('control',x,ignore.case = T)))]} 
condition.com <- data.matrix(condition.com)

dds <- DESeqDataSetFromMatrix(countData = final, colData = expdesign, design = ~condition)
dds <- DESeq(dds)
n <- ncol(condition.com) 
dir.create('/scratch/beoungle/probing_rnaseq/c42b/output/c42b/DEresults_DEseq2') 

for (i in 1:n) {
res <- results(dds, contrast = c('condition', condition.com[1,i], condition.com[2,i]))
res <- res[order(res$padj),]
upregulated <- subset(res, padj < 0.05 & (log2FoldChange > 0.5849625))
downregulated <- subset(res, padj < 0.05 & (log2FoldChange < -0.5849625)) 
if(nrow(upregulated) != 0) {write.csv(upregulated, file = paste0('/scratch/beoungle/probing_rnaseq/c42b/output/c42b/DEresults_DEseq2/', condition.com[1,i],'.VS.', condition.com[2,i],'.upregulated.csv'), quote = F, row.names = T)}
if(nrow(downregulated) != 0) {write.csv(downregulated, file = paste0('/scratch/beoungle/probing_rnaseq/c42b/output/c42b/DEresults_DEseq2/', condition.com[1,i],'.VS.', condition.com[2,i],'.downregulated.csv'), quote = F, row.names = T)}
if(nrow(res) != 0) {write.csv(res, file = paste0('/scratch/beoungle/probing_rnaseq/c42b/output/c42b/DEresults_DEseq2/', condition.com[1,i],'.VS.', condition.com[2,i],'.ALL.csv'), quote = F, row.names = T)}
}