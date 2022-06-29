library(dplyr)
library(ggplot2)




reading_deseq <- read.table("C42B_mean_fpkm_greater_than_0.5.tsv", skip = 1, row.names = NULL, sep = "\t")

reading_gencode <- read.table("Promoter_in_loop_with_promoter.chr.bed")

colnames(reading_deseq)[1] <- "gene_id"
colnames(reading_gencode)[4] <- "gene_id"

promoter_merged_deseq <- merge(reading_deseq, reading_gencode, by = "gene_id")


reading_gencode <- read.table("Promoter_in_loop_with_enhancer.chr.bed")

colnames(reading_deseq)[1] <- "gene_id"
colnames(reading_gencode)[4] <- "gene_id"

enhancer_merged_deseq <- merge(reading_deseq, reading_gencode, by = "gene_id")

reading_gencode <- read.table("Promoter_in_loop_with_insulator.chr.bed")

colnames(reading_gencode)[4] <- "gene_id"

insulator_merged_deseq <- merge(reading_deseq, reading_gencode, by = "gene_id")

reading_gencode <- read.table("Promoter_in_loop_with_ndr.chr.bed")

colnames(reading_deseq)[1] <- "gene_id"
colnames(reading_gencode)[4] <- "gene_id"

ndr_merged_deseq <- merge(reading_deseq, reading_gencode, by = "gene_id")


reading_gencode <- read.table("Promoter_in_loop_with_H3K27me3.chr.bed")

colnames(reading_deseq)[1] <- "gene_id"
colnames(reading_gencode)[4] <- "gene_id"

h3k27me3_merged_deseq <- merge(reading_deseq, reading_gencode, by = "gene_id")

reading_gencode <- read.table("Promoter_in_loop_with_h3k9me3.chr.bed")

colnames(reading_deseq)[1] <- "gene_id"
colnames(reading_gencode)[4] <- "gene_id"

h3k9me3_merged_deseq <- merge(reading_deseq, reading_gencode, by = "gene_id")


reading_gencode <- read.table("Promoter_in_loop_with_none.chr.bed")

colnames(reading_deseq)[1] <- "gene_id"
colnames(reading_gencode)[4] <- "gene_id"

none_merged_deseq <- merge(reading_deseq, reading_gencode, by = "gene_id")

promoter_merged_deseq['Loop_Category'] <- 'Promoter-Promoter'
enhancer_merged_deseq['Loop_Category'] <- "Promoter-Enhancer"
insulator_merged_deseq['Loop_Category'] <- "Promoter-Insulator"
ndr_merged_deseq['Loop_Category'] <- "Promoter-NDR"
h3k27me3_merged_deseq['Loop_Category'] <- "Promoter-H3K27me3"
h3k9me3_merged_deseq['Loop_Category'] <- "Promoter-H3K9me3"
none_merged_deseq['Loop_Category'] <- "Promoter-None"


merging_deseq <- rbind(promoter_merged_deseq, enhancer_merged_deseq, insulator_merged_deseq, ndr_merged_deseq, h3k27me3_merged_deseq, h3k9me3_merged_deseq, none_merged_deseq)

jpeg("Promoter_itself_by_loop_category.jpg", width = 1600, height = 1600)
ggplot(data = merging_deseq, aes(y=V23, x=Loop_Category)) +geom_violin( fill = "grey") + ylim(0,25) + ylab("FPKM") + stat_summary(fun = "mean",geom = "point",color = "red", size = 5) + stat_summary(fun = "median", geom = "point",color = "blue", size = 5)+ xlab("Promoter Categories") + theme(axis.text =  element_text(size = 10)) + scale_fill_grey() + theme_bw() + theme(  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
dev.off()




#all_promoter_p <- c(t.test(ndr_merged_deseq$V23, enhancer_merged_deseq$V23)$p.value, t.test(ndr_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(ndr_merged_deseq$V23, insulator_merged_deseq$V23)$p.value, t.test(ndr_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(ndr_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, insulator_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(h3k9me3_merged_deseq$V23, insulator_merged_deseq$V23)$p.value, t.test(h3k9me3_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(h3k9me3_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(insulator_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(insulator_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(h3k27me3_merged_deseq$V23, none_merged_deseq$V23)$p.value)

#p.adjust(all_promoter_p, method ="fdr", n = length(all_promoter_p))

#all_promoter_p <- c(t.test(enhancer_merged_deseq$V23, insulator_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, ndr_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(insulator_merged_deseq$V23, ndr_merged_deseq$V23)$p.value, t.test(insulator_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(insulator_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(insulator_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(ndr_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(ndr_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(ndr_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(h3k27me3_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(h3k27me3_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(h3k9me3_merged_deseq$V23, none_merged_deseq$V23)$p.value)

all_promoter_p <- c(t.test(promoter_merged_deseq$V23, enhancer_merged_deseq$V23)$p.value,t.test(promoter_merged_deseq$V23, insulator_merged_deseq$V23)$p.value, t.test(promoter_merged_deseq$V23, ndr_merged_deseq$V23)$p.value, t.test(promoter_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(promoter_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(promoter_merged_deseq$V23, none_merged_deseq$V23)$p.value,t.test(enhancer_merged_deseq$V23, insulator_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, ndr_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(enhancer_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(insulator_merged_deseq$V23, ndr_merged_deseq$V23)$p.value, t.test(insulator_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(insulator_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(insulator_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(ndr_merged_deseq$V23, h3k27me3_merged_deseq$V23)$p.value, t.test(ndr_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(ndr_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(h3k27me3_merged_deseq$V23, h3k9me3_merged_deseq$V23)$p.value, t.test(h3k27me3_merged_deseq$V23, none_merged_deseq$V23)$p.value, t.test(h3k9me3_merged_deseq$V23, none_merged_deseq$V23)$p.value)


p.adjust(all_promoter_p, method ="fdr", n = length(all_promoter_p))

