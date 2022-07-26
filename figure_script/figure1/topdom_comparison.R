library(dplyr)
library(ggplot2)
library(plyr)
library(fuzzyjoin)

method_hic <- read.table("./input/HiC-50000-topdom.true.bed")
method_micro <- read.table("./input/1_Billion-50000-topdom.ture.bed")

method_micro$V1 <- gsub("chr","",method_micro$V1)

method_hic$V1 <- gsub("chr","",method_hic$V1)

method_hic <- method_hic %>% mutate( MicroC = 0,Capture_Micro_C = 1)
method_micro <- method_micro %>% mutate (MicroC = 1,Capture_Micro_C = 0)


method_hic$loopname <- paste("method_hic",rownames(method_hic),sep="_")
method_micro$loopname <-paste("method_micro", rownames(method_micro),sep="_")


method_hic <- method_hic %>% filter(V4 == "domain")
method_micro <- method_micro %>% filter(V4 == "domain")

x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")



micro_2_hic <- data.frame()
for (i in x){
  chr_template <- method_micro %>% filter(V1 == i)
  chr_reference <- method_hic %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 100000, by = c("V2","V3"))
  micro_2_hic <-rbind(micro_2_hic,some_count)
}
method_micro$Capture_Micro_C[method_micro$loopname %in% micro_2_hic$loopname] <- 1

nrow(hic_micro_2)
