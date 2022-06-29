
library(HiCcompare)
library(TopDom)
library(readr)

#Run TopDom

chr_list <- c("chr1","chr2","chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

args <- commandArgs(trailingOnly=TRUE)

prefix_name <- args[1]
resolution <- args[2]

options(scipen = 99) #Force to not use scietific notation

for ( i in chr_list) {
matrix_list <- paste0(prefix_name,"_",i,"_",resolution,"bp.txt")
matrix_sparse <- read.table(matrix_list, header = FALSE)
mat <- sparse2full(matrix_sparse)
bed <- data.frame( chr= i, start = colnames(mat),end = as.numeric(colnames(mat)) + 50000)
mat <- cbind(bed, mat)
matrix_tsv <- paste0(prefix_name,"_",i,"_",resolution,"bp.matrix",sep="")
write_tsv(mat, file = matrix_tsv, col_names = FALSE)
topdom_out <- paste0(prefix_name,"_",i,"_",resolution,"bp.topdom",sep="")
TopDom(data = matrix_tsv, window.size = 5,outFile=topdom_out)
}



