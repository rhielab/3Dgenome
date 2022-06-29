library(plyr)
library(dplyr)
library(fuzzyjoin)
library(ggplot2)
library(ggpubr)
library(UpSetR)

#This script is similar to loop category script, but filter loops based on the loop information they have, and save them based on the loop categories so we can compare the loop information.

#I may be able to combine some of these scripts for paper purpose?

c42b_micro_2 <- read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/2_Billion_Mustache-5kb-all-loop.bedpe")
c42b_chicago <-  read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/PC-UNI1945.5kb.bedpe")
c42b_hic <-  read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/HiC_Mustache-5kb-all-loop.bedpe")
c42b_micro_1 <- read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/1_Billion_Mustache-5kb-all-loop.bedpe")
c42b_micro_3 <- read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/3_billion_low_Mustache-5kb-all-loop.bedpe")

c42b_micro_2$loopname <- paste("c42b_2_billion",rownames(c42b_micro_2),sep = "_")
c42b_micro_1$loopname <- paste("c42b_1_billion",rownames(c42b_micro_1),sep = "_")
c42b_chicago$loopname <- paste("c42b_chicago",rownames(c42b_chicago),sep = "_")
c42b_hic$loopname <- paste("c42b_hic",rownames(c42b_hic),sep = "_")
c42b_micro_3$loopname <- paste("c42b_3_billion",rownames(c42b_micro_3),sep = "_")


h3k27ac <- read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/C42B_H3K27ac_no_tss.bed")
h3k4me3 <- read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/C42B_TSS_FPKM_greater_than_0.5.tsv")
ctcf <- read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/C42B_CTCF_no_tss_no_enhancer.bed")
ndr <- read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/C42B_NDR_without_features.bed")
h3k27me3 <- read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/H3K27me3_not_in_enhancer_insulator_promoter_ndr_region.bed")
h3k9me3 <- read.table("/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/H3K9me3_not_in_enhancer_insulator_promoter_ndr_h3k27me3_region.bed")

ctcf$V1 <- gsub("chr","",ctcf$V1) 
h3k27ac$V1 <- gsub("chr","",h3k27ac$V1)
h3k4me3$V1 <- gsub("chr","",h3k4me3$V1)
ndr$V1 <- gsub("chr","",ndr$V1)
h3k27me3$V1 <- gsub("chr","",h3k27me3$V1)
h3k9me3$V1 <- gsub("chr","",h3k9me3$V1)

c42b_micro_2$V1 <- gsub("chr","",c42b_micro_2$V1)
c42b_micro_2$V4 <- gsub("chr","",c42b_micro_2$V4)

c42b_micro_1$V1 <- gsub("chr","",c42b_micro_1$V1)
c42b_micro_1$V4 <- gsub("chr","",c42b_micro_1$V4)

c42b_micro_3$V1 <- gsub("chr","",c42b_micro_3$V1)
c42b_micro_3$V4 <- gsub("chr","",c42b_micro_3$V4)

c42b_hic$V1 <- gsub("chr","",c42b_hic$V1)
c42b_hic$V4 <- gsub("chr","",c42b_hic$V4)


head(h3k4me3)

c42b_micro_2 <- c42b_micro_2 %>% mutate(Enhancer_1 = 0, Enhancer_2 = 0, Promoter_1 = 0, Promoter_2 = 0, Insulator_1 = 0, Insulator_2 = 0, NDR_1 = 0, NDR_2 = 0, H3K27me3_1 = 0, H3K27me3_2 = 0, H3K9me3_1 = 0, H3K9me3_2 = 0, None_1 = 1, None_2 = 1)
c42b_hic <- c42b_hic %>% mutate(Enhancer_1 = 0, Enhancer_2 = 0, Promoter_1 = 0, Promoter_2 = 0, Insulator_1 = 0, Insulator_2 = 0, NDR_1 = 0, NDR_2 = 0, H3K27me3_1 = 0, H3K27me3_2 = 0, H3K9me3_1 = 0, H3K9me3_2 = 0, None_1 = 1, None_2 = 1)
c42b_micro_1 <- c42b_micro_1 %>% mutate(Enhancer_1 = 0, Enhancer_2 = 0, Promoter_1 = 0, Promoter_2 = 0, Insulator_1 = 0, Insulator_2 = 0, NDR_1 = 0, NDR_2 = 0, H3K27me3_1 = 0, H3K27me3_2 = 0, H3K9me3_1 = 0, H3K9me3_2 = 0, None_1 = 1, None_2 = 1)
c42b_chicago <- c42b_chicago %>% mutate(Enhancer_1 = 0, Enhancer_2 = 0, Promoter_1 = 0, Promoter_2 = 0, Insulator_1 = 0, Insulator_2 = 0, NDR_1 = 0, NDR_2 = 0, H3K27me3_1 = 0, H3K27me3_2 = 0, H3K9me3_1 = 0, H3K9me3_2 = 0, None_1 = 1, None_2 = 1)
c42b_micro_3 <- c42b_micro_3 %>% mutate(Enhancer_1 = 0, Enhancer_2 = 0, Promoter_1 = 0, Promoter_2 = 0, Insulator_1 = 0, Insulator_2 = 0, NDR_1 = 0, NDR_2 = 0, H3K27me3_1 = 0, H3K27me3_2 = 0, H3K9me3_1 = 0, H3K9me3_2 = 0, None_1 = 1, None_2 = 1)

x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")

#Micro C 2 Billion_5kb


c42b_micro_2_anchor1 <- c42b_micro_2 %>% select(V1,V2,V3,loopname)
c42b_micro_2_anchor2 <- c42b_micro_2 %>% select(V4,V5,V6,loopname)

colnames(c42b_micro_2_anchor2)[1] <- "V1"
colnames(c42b_micro_2_anchor2)[2] <- "V2"
colnames(c42b_micro_2_anchor2)[3] <- "V3"


c42b_micro_2_1_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_2_anchor1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_1_h3k4me3 <-rbind(c42b_micro_2_1_h3k4me3,some_count)
}

c42b_micro_2_2_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_2_anchor2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_2_h3k4me3 <-rbind(c42b_micro_2_2_h3k4me3,some_count)
}


c42b_micro_2$Promoter_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_h3k4me3$loopname] <- 1
c42b_micro_2$Promoter_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_h3k4me3$loopname] <- 1
c42b_micro_2$None_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_h3k4me3$loopname] <- 0
c42b_micro_2$None_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_h3k4me3$loopname] <- 0

head(c42b_micro_2)
head(c42b_micro_2_1_h3k4me3)

non_c42b_micro_2_1_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_2_anchor1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_2_1_h3k4me3 <-rbind(non_c42b_micro_2_1_h3k4me3,some_count)
}

non_c42b_micro_2_2_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_2_anchor2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_2_2_h3k4me3 <-rbind(non_c42b_micro_2_2_h3k4me3,some_count)
}


c42b_micro_2_1_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_1_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_1_h3k27ac <-rbind(c42b_micro_2_1_h3k27ac,some_count)
}

c42b_micro_2_2_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_2_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_2_h3k27ac <-rbind(c42b_micro_2_2_h3k27ac,some_count)
}

c42b_micro_2$Enhancer_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_h3k27ac$loopname] <- 1
c42b_micro_2$Enhancer_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_h3k27ac$loopname] <- 1
c42b_micro_2$None_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_h3k27ac$loopname] <- 0
c42b_micro_2$None_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_h3k27ac$loopname] <- 0


non_c42b_micro_2_1_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_1_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_2_1_h3k27ac <-rbind(non_c42b_micro_2_1_h3k27ac,some_count)
}

non_c42b_micro_2_2_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_2_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_2_2_h3k27ac <-rbind(non_c42b_micro_2_2_h3k27ac,some_count)
}


c42b_micro_2_1_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_1_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_1_ctcf <-rbind(c42b_micro_2_1_ctcf,some_count)
}


c42b_micro_2_2_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_2_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_2_ctcf <-rbind(c42b_micro_2_2_ctcf,some_count)
}

c42b_micro_2$Insulator_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_ctcf$loopname] <- 1
c42b_micro_2$Insulator_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_ctcf$loopname] <- 1
c42b_micro_2$None_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_ctcf$loopname] <- 0
c42b_micro_2$None_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_ctcf$loopname] <- 0

non_c42b_micro_2_1_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_1_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_2_1_ctcf <-rbind(non_c42b_micro_2_1_ctcf,some_count)
}


non_c42b_micro_2_2_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_2_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_2_2_ctcf <-rbind(non_c42b_micro_2_2_ctcf,some_count)
}


c42b_micro_2_1_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_1_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_1_ndr <-rbind(c42b_micro_2_1_ndr,some_count)
}

c42b_micro_2_2_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_2_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_2_ndr <-rbind(c42b_micro_2_2_ndr,some_count)
}


c42b_micro_2$NDR_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_ndr$loopname] <- 1
c42b_micro_2$NDR_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_ndr$loopname] <- 1
c42b_micro_2$None_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_ndr$loopname] <- 0
c42b_micro_2$None_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_ndr$loopname] <- 0

non_c42b_micro_2_1_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_1_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_2_1_ndr <-rbind(non_c42b_micro_2_1_ndr,some_count)
}

non_c42b_micro_2_2_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_2_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_2_2_ndr <-rbind(non_c42b_micro_2_2_ndr,some_count)
}

c42b_micro_2_1_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_1_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_1_h3k27me3 <-rbind(c42b_micro_2_1_h3k27me3,some_count)
}

c42b_micro_2_2_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_2_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_2_h3k27me3 <-rbind(c42b_micro_2_2_h3k27me3,some_count)
}


c42b_micro_2$H3K27me3_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_h3k27me3$loopname] <- 1
c42b_micro_2$H3K27me3_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_h3k27me3$loopname] <- 1
c42b_micro_2$None_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_h3k27me3$loopname] <- 0
c42b_micro_2$None_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_h3k27me3$loopname] <- 0


non_c42b_micro_2_1_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_1_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_2_1_h3k27me3 <-rbind(non_c42b_micro_2_1_h3k27me3,some_count)
}

non_c42b_micro_2_2_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_2_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_2_2_h3k27me3 <-rbind(non_c42b_micro_2_2_h3k27me3,some_count)
}

c42b_micro_2_1_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_1_h3k27me3 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_1_h3k9me3 <-rbind(c42b_micro_2_1_h3k9me3,some_count)
}

c42b_micro_2_2_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_2_2_h3k27me3 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_2_2_h3k9me3 <-rbind(c42b_micro_2_2_h3k9me3,some_count)
}


c42b_micro_2$H3K9me3_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_h3k9me3$loopname] <- 1
c42b_micro_2$H3K9me3_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_h3k9me3$loopname] <- 1
c42b_micro_2$None_1[c42b_micro_2$loopname  %in% c42b_micro_2_1_h3k9me3$loopname] <- 0
c42b_micro_2$None_2[c42b_micro_2$loopname  %in% c42b_micro_2_2_h3k9me3$loopname] <- 0


c42b_micro_1_anchor1 <- c42b_micro_1 %>% select(V1,V2,V3,loopname)
c42b_micro_1_anchor2 <- c42b_micro_1 %>% select(V4,V5,V6,loopname)

colnames(c42b_micro_1_anchor2)[1] <- "V1"
colnames(c42b_micro_1_anchor2)[2] <- "V2"
colnames(c42b_micro_1_anchor2)[3] <- "V3"


c42b_micro_1_1_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_1_anchor1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_1_h3k4me3 <-rbind(c42b_micro_1_1_h3k4me3,some_count)
}

c42b_micro_1_2_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_1_anchor2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_2_h3k4me3 <-rbind(c42b_micro_1_2_h3k4me3,some_count)
}


c42b_micro_1$Promoter_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_h3k4me3$loopname] <- 1
c42b_micro_1$Promoter_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_h3k4me3$loopname] <- 1
c42b_micro_1$None_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_h3k4me3$loopname] <- 0
c42b_micro_1$None_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_h3k4me3$loopname] <- 0

non_c42b_micro_1_1_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_1_anchor1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_1_1_h3k4me3 <-rbind(non_c42b_micro_1_1_h3k4me3,some_count)
}

non_c42b_micro_1_2_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_1_anchor2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_1_2_h3k4me3 <-rbind(non_c42b_micro_1_2_h3k4me3,some_count)
}


c42b_micro_1_1_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_1_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_1_h3k27ac <-rbind(c42b_micro_1_1_h3k27ac,some_count)
}

c42b_micro_1_2_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_2_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_2_h3k27ac <-rbind(c42b_micro_1_2_h3k27ac,some_count)
}

c42b_micro_1$Enhancer_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_h3k27ac$loopname] <- 1
c42b_micro_1$Enhancer_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_h3k27ac$loopname] <- 1
c42b_micro_1$None_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_h3k27ac$loopname] <- 0
c42b_micro_1$None_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_h3k27ac$loopname] <- 0


non_c42b_micro_1_1_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_1_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_1_1_h3k27ac <-rbind(non_c42b_micro_1_1_h3k27ac,some_count)
}

non_c42b_micro_1_2_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_2_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_1_2_h3k27ac <-rbind(non_c42b_micro_1_2_h3k27ac,some_count)
}


c42b_micro_1_1_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_1_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_1_ctcf <-rbind(c42b_micro_1_1_ctcf,some_count)
}


c42b_micro_1_2_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_2_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_2_ctcf <-rbind(c42b_micro_1_2_ctcf,some_count)
}

c42b_micro_1$Insulator_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_ctcf$loopname] <- 1
c42b_micro_1$Insulator_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_ctcf$loopname] <- 1
c42b_micro_1$None_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_ctcf$loopname] <- 0
c42b_micro_1$None_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_ctcf$loopname] <- 0

non_c42b_micro_1_1_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_1_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_1_1_ctcf <-rbind(non_c42b_micro_1_1_ctcf,some_count)
}


non_c42b_micro_1_2_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_2_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_1_2_ctcf <-rbind(non_c42b_micro_1_2_ctcf,some_count)
}


c42b_micro_1_1_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_1_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_1_ndr <-rbind(c42b_micro_1_1_ndr,some_count)
}

c42b_micro_1_2_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_2_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_2_ndr <-rbind(c42b_micro_1_2_ndr,some_count)
}


c42b_micro_1$NDR_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_ndr$loopname] <- 1
c42b_micro_1$NDR_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_ndr$loopname] <- 1
c42b_micro_1$None_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_ndr$loopname] <- 0
c42b_micro_1$None_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_ndr$loopname] <- 0

non_c42b_micro_1_1_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_1_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_1_1_ndr <-rbind(non_c42b_micro_1_1_ndr,some_count)
}

non_c42b_micro_1_2_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_2_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_1_2_ndr <-rbind(non_c42b_micro_1_2_ndr,some_count)
}

c42b_micro_1_1_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_1_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_1_h3k27me3 <-rbind(c42b_micro_1_1_h3k27me3,some_count)
}

c42b_micro_1_2_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_2_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_2_h3k27me3 <-rbind(c42b_micro_1_2_h3k27me3,some_count)
}


c42b_micro_1$H3K27me3_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_h3k27me3$loopname] <- 1
c42b_micro_1$H3K27me3_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_h3k27me3$loopname] <- 1
c42b_micro_1$None_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_h3k27me3$loopname] <- 0
c42b_micro_1$None_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_h3k27me3$loopname] <- 0


non_c42b_micro_1_1_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_1_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_1_1_h3k27me3 <-rbind(non_c42b_micro_1_1_h3k27me3,some_count)
}

non_c42b_micro_1_2_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_2_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_1_2_h3k27me3 <-rbind(non_c42b_micro_1_2_h3k27me3,some_count)
}

c42b_micro_1_1_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_1_h3k27me3 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_1_h3k9me3 <-rbind(c42b_micro_1_1_h3k9me3,some_count)
}

c42b_micro_1_2_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_1_2_h3k27me3 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_1_2_h3k9me3 <-rbind(c42b_micro_1_2_h3k9me3,some_count)
}


c42b_micro_1$H3K9me3_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_h3k9me3$loopname] <- 1
c42b_micro_1$H3K9me3_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_h3k9me3$loopname] <- 1
c42b_micro_1$None_1[c42b_micro_1$loopname  %in% c42b_micro_1_1_h3k9me3$loopname] <- 0
c42b_micro_1$None_2[c42b_micro_1$loopname  %in% c42b_micro_1_2_h3k9me3$loopname] <- 0


c42b_micro_3_anchor1 <- c42b_micro_3 %>% select(V1,V2,V3,loopname)
c42b_micro_3_anchor2 <- c42b_micro_3 %>% select(V4,V5,V6,loopname)

colnames(c42b_micro_3_anchor2)[1] <- "V1"
colnames(c42b_micro_3_anchor2)[2] <- "V2"
colnames(c42b_micro_3_anchor2)[3] <- "V3"


c42b_micro_3_1_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_anchor1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_1_h3k4me3 <-rbind(c42b_micro_3_1_h3k4me3,some_count)
}

c42b_micro_3_2_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_anchor2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_2_h3k4me3 <-rbind(c42b_micro_3_2_h3k4me3,some_count)
}


c42b_micro_3$Promoter_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_h3k4me3$loopname] <- 1
c42b_micro_3$Promoter_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_h3k4me3$loopname] <- 1
c42b_micro_3$None_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_h3k4me3$loopname] <- 0
c42b_micro_3$None_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_h3k4me3$loopname] <- 0

non_c42b_micro_3_1_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_anchor1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_3_1_h3k4me3 <-rbind(non_c42b_micro_3_1_h3k4me3,some_count)
}

non_c42b_micro_3_2_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_anchor2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_3_2_h3k4me3 <-rbind(non_c42b_micro_3_2_h3k4me3,some_count)
}


c42b_micro_3_1_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_1_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_1_h3k27ac <-rbind(c42b_micro_3_1_h3k27ac,some_count)
}

c42b_micro_3_2_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_2_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_2_h3k27ac <-rbind(c42b_micro_3_2_h3k27ac,some_count)
}

c42b_micro_3$Enhancer_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_h3k27ac$loopname] <- 1
c42b_micro_3$Enhancer_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_h3k27ac$loopname] <- 1
c42b_micro_3$None_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_h3k27ac$loopname] <- 0
c42b_micro_3$None_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_h3k27ac$loopname] <- 0


non_c42b_micro_3_1_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_1_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_3_1_h3k27ac <-rbind(non_c42b_micro_3_1_h3k27ac,some_count)
}

non_c42b_micro_3_2_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_2_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_3_2_h3k27ac <-rbind(non_c42b_micro_3_2_h3k27ac,some_count)
}


c42b_micro_3_1_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_1_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_1_ctcf <-rbind(c42b_micro_3_1_ctcf,some_count)
}


c42b_micro_3_2_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_2_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_2_ctcf <-rbind(c42b_micro_3_2_ctcf,some_count)
}

c42b_micro_3$Insulator_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_ctcf$loopname] <- 1
c42b_micro_3$Insulator_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_ctcf$loopname] <- 1
c42b_micro_3$None_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_ctcf$loopname] <- 0
c42b_micro_3$None_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_ctcf$loopname] <- 0

non_c42b_micro_3_1_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_1_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_3_1_ctcf <-rbind(non_c42b_micro_3_1_ctcf,some_count)
}


non_c42b_micro_3_2_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_2_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_3_2_ctcf <-rbind(non_c42b_micro_3_2_ctcf,some_count)
}


c42b_micro_3_1_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_1_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_1_ndr <-rbind(c42b_micro_3_1_ndr,some_count)
}

c42b_micro_3_2_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_2_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_2_ndr <-rbind(c42b_micro_3_2_ndr,some_count)
}


c42b_micro_3$NDR_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_ndr$loopname] <- 1
c42b_micro_3$NDR_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_ndr$loopname] <- 1
c42b_micro_3$None_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_ndr$loopname] <- 0
c42b_micro_3$None_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_ndr$loopname] <- 0

non_c42b_micro_3_1_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_1_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_3_1_ndr <-rbind(non_c42b_micro_3_1_ndr,some_count)
}

non_c42b_micro_3_2_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_2_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_3_2_ndr <-rbind(non_c42b_micro_3_2_ndr,some_count)
}

c42b_micro_3_1_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_1_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_1_h3k27me3 <-rbind(c42b_micro_3_1_h3k27me3,some_count)
}

c42b_micro_3_2_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_2_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_2_h3k27me3 <-rbind(c42b_micro_3_2_h3k27me3,some_count)
}


c42b_micro_3$H3K27me3_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_h3k27me3$loopname] <- 1
c42b_micro_3$H3K27me3_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_h3k27me3$loopname] <- 1
c42b_micro_3$None_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_h3k27me3$loopname] <- 0
c42b_micro_3$None_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_h3k27me3$loopname] <- 0


non_c42b_micro_3_1_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_1_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_3_1_h3k27me3 <-rbind(non_c42b_micro_3_1_h3k27me3,some_count)
}

non_c42b_micro_3_2_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_2_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_micro_3_2_h3k27me3 <-rbind(non_c42b_micro_3_2_h3k27me3,some_count)
}

c42b_micro_3_1_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_1_h3k27me3 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_1_h3k9me3 <-rbind(c42b_micro_3_1_h3k9me3,some_count)
}

c42b_micro_3_2_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_micro_3_2_h3k27me3 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_2_h3k9me3 <-rbind(c42b_micro_3_2_h3k9me3,some_count)
}


c42b_micro_3$H3K9me3_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_h3k9me3$loopname] <- 1
c42b_micro_3$H3K9me3_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_h3k9me3$loopname] <- 1
c42b_micro_3$None_1[c42b_micro_3$loopname  %in% c42b_micro_3_1_h3k9me3$loopname] <- 0
c42b_micro_3$None_2[c42b_micro_3$loopname  %in% c42b_micro_3_2_h3k9me3$loopname] <- 0


c42b_hic_anchor1 <- c42b_hic %>% select(V1,V2,V3,loopname)
c42b_hic_anchor2 <- c42b_hic %>% select(V4,V5,V6,loopname)

colnames(c42b_hic_anchor2)[1] <- "V1"
colnames(c42b_hic_anchor2)[2] <- "V2"
colnames(c42b_hic_anchor2)[3] <- "V3"


c42b_hic_1_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_hic_anchor1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_1_h3k4me3 <-rbind(c42b_hic_1_h3k4me3,some_count)
}

c42b_hic_2_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_hic_anchor2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_2_h3k4me3 <-rbind(c42b_hic_2_h3k4me3,some_count)
}


c42b_hic$Promoter_1[c42b_hic$loopname  %in% c42b_hic_1_h3k4me3$loopname] <- 1
c42b_hic$Promoter_2[c42b_hic$loopname  %in% c42b_hic_2_h3k4me3$loopname] <- 1
c42b_hic$None_1[c42b_hic$loopname  %in% c42b_hic_1_h3k4me3$loopname] <- 0
c42b_hic$None_2[c42b_hic$loopname  %in% c42b_hic_2_h3k4me3$loopname] <- 0

non_c42b_hic_1_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_hic_anchor1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_hic_1_h3k4me3 <-rbind(non_c42b_hic_1_h3k4me3,some_count)
}

non_c42b_hic_2_h3k4me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_hic_anchor2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_hic_2_h3k4me3 <-rbind(non_c42b_hic_2_h3k4me3,some_count)
}


c42b_hic_1_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_1_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_1_h3k27ac <-rbind(c42b_hic_1_h3k27ac,some_count)
}

c42b_hic_2_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_2_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_2_h3k27ac <-rbind(c42b_hic_2_h3k27ac,some_count)
}

c42b_hic$Enhancer_1[c42b_hic$loopname  %in% c42b_hic_1_h3k27ac$loopname] <- 1
c42b_hic$Enhancer_2[c42b_hic$loopname  %in% c42b_hic_2_h3k27ac$loopname] <- 1
c42b_hic$None_1[c42b_hic$loopname  %in% c42b_hic_1_h3k27ac$loopname] <- 0
c42b_hic$None_2[c42b_hic$loopname  %in% c42b_hic_2_h3k27ac$loopname] <- 0


non_c42b_hic_1_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_1_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_hic_1_h3k27ac <-rbind(non_c42b_hic_1_h3k27ac,some_count)
}

non_c42b_hic_2_h3k27ac <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_2_h3k4me3 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_hic_2_h3k27ac <-rbind(non_c42b_hic_2_h3k27ac,some_count)
}


c42b_hic_1_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_1_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_1_ctcf <-rbind(c42b_hic_1_ctcf,some_count)
}


c42b_hic_2_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_2_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_2_ctcf <-rbind(c42b_hic_2_ctcf,some_count)
}

c42b_hic$Insulator_1[c42b_hic$loopname  %in% c42b_hic_1_ctcf$loopname] <- 1
c42b_hic$Insulator_2[c42b_hic$loopname  %in% c42b_hic_2_ctcf$loopname] <- 1
c42b_hic$None_1[c42b_hic$loopname  %in% c42b_hic_1_ctcf$loopname] <- 0
c42b_hic$None_2[c42b_hic$loopname  %in% c42b_hic_2_ctcf$loopname] <- 0

non_c42b_hic_1_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_1_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_hic_1_ctcf <-rbind(non_c42b_hic_1_ctcf,some_count)
}


non_c42b_hic_2_ctcf <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_2_h3k27ac %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_hic_2_ctcf <-rbind(non_c42b_hic_2_ctcf,some_count)
}


c42b_hic_1_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_1_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_1_ndr <-rbind(c42b_hic_1_ndr,some_count)
}

c42b_hic_2_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_2_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_2_ndr <-rbind(c42b_hic_2_ndr,some_count)
}


c42b_hic$NDR_1[c42b_hic$loopname  %in% c42b_hic_1_ndr$loopname] <- 1
c42b_hic$NDR_2[c42b_hic$loopname  %in% c42b_hic_2_ndr$loopname] <- 1
c42b_hic$None_1[c42b_hic$loopname  %in% c42b_hic_1_ndr$loopname] <- 0
c42b_hic$None_2[c42b_hic$loopname  %in% c42b_hic_2_ndr$loopname] <- 0

non_c42b_hic_1_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_1_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_hic_1_ndr <-rbind(non_c42b_hic_1_ndr,some_count)
}

non_c42b_hic_2_ndr <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_2_ctcf %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_hic_2_ndr <-rbind(non_c42b_hic_2_ndr,some_count)
}

c42b_hic_1_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_1_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_1_h3k27me3 <-rbind(c42b_hic_1_h3k27me3,some_count)
}

c42b_hic_2_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_2_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_2_h3k27me3 <-rbind(c42b_hic_2_h3k27me3,some_count)
}


c42b_hic$H3K27me3_1[c42b_hic$loopname  %in% c42b_hic_1_h3k27me3$loopname] <- 1
c42b_hic$H3K27me3_2[c42b_hic$loopname  %in% c42b_hic_2_h3k27me3$loopname] <- 1
c42b_hic$None_1[c42b_hic$loopname  %in% c42b_hic_1_h3k27me3$loopname] <- 0
c42b_hic$None_2[c42b_hic$loopname  %in% c42b_hic_2_h3k27me3$loopname] <- 0


non_c42b_hic_1_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_1_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_hic_1_h3k27me3 <-rbind(non_c42b_hic_1_h3k27me3,some_count)
}

non_c42b_hic_2_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_2_ndr %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_c42b_hic_2_h3k27me3 <-rbind(non_c42b_hic_2_h3k27me3,some_count)
}

c42b_hic_1_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_1_h3k27me3 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_1_h3k9me3 <-rbind(c42b_hic_1_h3k9me3,some_count)
}

c42b_hic_2_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- non_c42b_hic_2_h3k27me3 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  c42b_hic_2_h3k9me3 <-rbind(c42b_hic_2_h3k9me3,some_count)
}


c42b_hic$H3K9me3_1[c42b_hic$loopname  %in% c42b_hic_1_h3k9me3$loopname] <- 1
c42b_hic$H3K9me3_2[c42b_hic$loopname  %in% c42b_hic_2_h3k9me3$loopname] <- 1
c42b_hic$None_1[c42b_hic$loopname  %in% c42b_hic_1_h3k9me3$loopname] <- 0
c42b_hic$None_2[c42b_hic$loopname  %in% c42b_hic_2_h3k9me3$loopname] <- 0



nrow(filter(c42b_hic, Promoter_1 == 1, Promoter_2 ==1))
nrow(filter(c42b_hic, Promoter_1 == 1, Enhancer_2 ==1)) + nrow(filter(c42b_hic, Promoter_2 == 1, Enhancer_1 ==1))
nrow(filter(c42b_hic, Promoter_1 == 1, Insulator_2 ==1)) + nrow(filter(c42b_hic, Promoter_2 == 1, Insulator_1 ==1))
nrow(filter(c42b_hic, Promoter_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_hic, Promoter_2 == 1, NDR_1 ==1))
nrow(filter(c42b_hic, Promoter_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_hic, Promoter_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_hic, Promoter_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_hic, Promoter_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_hic, Promoter_1 == 1, None_2 ==1)) + nrow(filter(c42b_hic, Promoter_2 == 1, None_1 ==1))


nrow(filter(c42b_hic, Enhancer_1 == 1, Enhancer_2 ==1))
nrow(filter(c42b_hic, Enhancer_1 == 1, Insulator_2 ==1)) + nrow(filter(c42b_hic, Enhancer_2 == 1, Insulator_1 ==1))
nrow(filter(c42b_hic, Enhancer_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_hic, Enhancer_2 == 1, NDR_1 ==1))
nrow(filter(c42b_hic, Enhancer_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_hic, Enhancer_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_hic, Enhancer_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_hic, Enhancer_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_hic, Enhancer_1 == 1, None_2 ==1)) + nrow(filter(c42b_hic, Enhancer_2 == 1, None_1 ==1))


nrow(filter(c42b_hic, Insulator_1 == 1, Insulator_2 ==1))
nrow(filter(c42b_hic, Insulator_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_hic, Insulator_2 == 1, NDR_1 ==1))
nrow(filter(c42b_hic, Insulator_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_hic, Insulator_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_hic, Insulator_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_hic, Insulator_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_hic, Insulator_1 == 1, None_2 ==1)) + nrow(filter(c42b_hic, Insulator_2 == 1, None_1 ==1))


nrow(filter(c42b_hic, NDR_1 == 1, NDR_2 ==1))
nrow(filter(c42b_hic, NDR_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_hic, NDR_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_hic, NDR_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_hic, NDR_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_hic, NDR_1 == 1, None_2 ==1)) + nrow(filter(c42b_hic, NDR_2 == 1, None_1 ==1))


nrow(filter(c42b_hic, H3K27me3_1 == 1, H3K27me3_2 ==1))
nrow(filter(c42b_hic, H3K27me3_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_hic, H3K27me3_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_hic, H3K27me3_1 == 1, None_2 ==1)) + nrow(filter(c42b_hic, H3K27me3_2 == 1, None_1 ==1))


nrow(filter(c42b_hic, H3K9me3_1 == 1, H3K9me3_2 ==1))
nrow(filter(c42b_hic, H3K9me3_1 == 1, None_2 ==1)) + nrow(filter(c42b_hic, H3K9me3_2 == 1, None_1 ==1))

nrow(filter(c42b_hic, None_1 == 1, None_2 ==1))



nrow(filter(c42b_micro_1, Promoter_1 == 1, Promoter_2 ==1))
nrow(filter(c42b_micro_1, Promoter_1 == 1, Enhancer_2 ==1)) + nrow(filter(c42b_micro_1, Promoter_2 == 1, Enhancer_1 ==1))
nrow(filter(c42b_micro_1, Promoter_1 == 1, Insulator_2 ==1)) + nrow(filter(c42b_micro_1, Promoter_2 == 1, Insulator_1 ==1))
nrow(filter(c42b_micro_1, Promoter_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_micro_1, Promoter_2 == 1, NDR_1 ==1))
nrow(filter(c42b_micro_1, Promoter_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_1, Promoter_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_1, Promoter_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_1, Promoter_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_1, Promoter_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_1, Promoter_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_1, Enhancer_1 == 1, Enhancer_2 ==1))
nrow(filter(c42b_micro_1, Enhancer_1 == 1, Insulator_2 ==1)) + nrow(filter(c42b_micro_1, Enhancer_2 == 1, Insulator_1 ==1))
nrow(filter(c42b_micro_1, Enhancer_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_micro_1, Enhancer_2 == 1, NDR_1 ==1))
nrow(filter(c42b_micro_1, Enhancer_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_1, Enhancer_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_1, Enhancer_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_1, Enhancer_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_1, Enhancer_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_1, Enhancer_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_1, Insulator_1 == 1, Insulator_2 ==1))
nrow(filter(c42b_micro_1, Insulator_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_micro_1, Insulator_2 == 1, NDR_1 ==1))
nrow(filter(c42b_micro_1, Insulator_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_1, Insulator_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_1, Insulator_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_1, Insulator_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_1, Insulator_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_1, Insulator_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_1, NDR_1 == 1, NDR_2 ==1))
nrow(filter(c42b_micro_1, NDR_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_1, NDR_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_1, NDR_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_1, NDR_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_1, NDR_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_1, NDR_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_1, H3K27me3_1 == 1, H3K27me3_2 ==1))
nrow(filter(c42b_micro_1, H3K27me3_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_1, H3K27me3_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_1, H3K27me3_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_1, H3K27me3_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_1, H3K9me3_1 == 1, H3K9me3_2 ==1))
nrow(filter(c42b_micro_1, H3K9me3_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_1, H3K9me3_2 == 1, None_1 ==1))

nrow(filter(c42b_micro_1, None_1 == 1, None_2 ==1))



nrow(filter(c42b_micro_2, Promoter_1 == 1, Promoter_2 ==1))
nrow(filter(c42b_micro_2, Promoter_1 == 1, Enhancer_2 ==1)) + nrow(filter(c42b_micro_2, Promoter_2 == 1, Enhancer_1 ==1))
nrow(filter(c42b_micro_2, Promoter_1 == 1, Insulator_2 ==1)) + nrow(filter(c42b_micro_2, Promoter_2 == 1, Insulator_1 ==1))
nrow(filter(c42b_micro_2, Promoter_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_micro_2, Promoter_2 == 1, NDR_1 ==1))
nrow(filter(c42b_micro_2, Promoter_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_2, Promoter_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_2, Promoter_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_2, Promoter_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_2, Promoter_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_2, Promoter_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_2, Enhancer_1 == 1, Enhancer_2 ==1))
nrow(filter(c42b_micro_2, Enhancer_1 == 1, Insulator_2 ==1)) + nrow(filter(c42b_micro_2, Enhancer_2 == 1, Insulator_1 ==1))
nrow(filter(c42b_micro_2, Enhancer_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_micro_2, Enhancer_2 == 1, NDR_1 ==1))
nrow(filter(c42b_micro_2, Enhancer_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_2, Enhancer_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_2, Enhancer_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_2, Enhancer_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_2, Enhancer_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_2, Enhancer_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_2, Insulator_1 == 1, Insulator_2 ==1))
nrow(filter(c42b_micro_2, Insulator_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_micro_2, Insulator_2 == 1, NDR_1 ==1))
nrow(filter(c42b_micro_2, Insulator_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_2, Insulator_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_2, Insulator_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_2, Insulator_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_2, Insulator_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_2, Insulator_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_2, NDR_1 == 1, NDR_2 ==1))
nrow(filter(c42b_micro_2, NDR_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_2, NDR_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_2, NDR_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_2, NDR_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_2, NDR_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_2, NDR_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_2, H3K27me3_1 == 1, H3K27me3_2 ==1))
nrow(filter(c42b_micro_2, H3K27me3_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_2, H3K27me3_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_2, H3K27me3_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_2, H3K27me3_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_2, H3K9me3_1 == 1, H3K9me3_2 ==1))
nrow(filter(c42b_micro_2, H3K9me3_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_2, H3K9me3_2 == 1, None_1 ==1))

nrow(filter(c42b_micro_2, None_1 == 1, None_2 ==1))



nrow(filter(c42b_micro_3, Promoter_1 == 1, Promoter_2 ==1))
nrow(filter(c42b_micro_3, Promoter_1 == 1, Enhancer_2 ==1)) + nrow(filter(c42b_micro_3, Promoter_2 == 1, Enhancer_1 ==1))
nrow(filter(c42b_micro_3, Promoter_1 == 1, Insulator_2 ==1)) + nrow(filter(c42b_micro_3, Promoter_2 == 1, Insulator_1 ==1))
nrow(filter(c42b_micro_3, Promoter_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_micro_3, Promoter_2 == 1, NDR_1 ==1))
nrow(filter(c42b_micro_3, Promoter_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_3, Promoter_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_3, Promoter_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_3, Promoter_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_3, Promoter_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_3, Promoter_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_3, Enhancer_1 == 1, Enhancer_2 ==1))
nrow(filter(c42b_micro_3, Enhancer_1 == 1, Insulator_2 ==1)) + nrow(filter(c42b_micro_3, Enhancer_2 == 1, Insulator_1 ==1))
nrow(filter(c42b_micro_3, Enhancer_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_micro_3, Enhancer_2 == 1, NDR_1 ==1))
nrow(filter(c42b_micro_3, Enhancer_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_3, Enhancer_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_3, Enhancer_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_3, Enhancer_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_3, Enhancer_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_3, Enhancer_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_3, Insulator_1 == 1, Insulator_2 ==1))
nrow(filter(c42b_micro_3, Insulator_1 == 1, NDR_2 ==1)) + nrow(filter(c42b_micro_3, Insulator_2 == 1, NDR_1 ==1))
nrow(filter(c42b_micro_3, Insulator_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_3, Insulator_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_3, Insulator_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_3, Insulator_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_3, Insulator_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_3, Insulator_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_3, NDR_1 == 1, NDR_2 ==1))
nrow(filter(c42b_micro_3, NDR_1 == 1, H3K27me3_2 ==1)) + nrow(filter(c42b_micro_3, NDR_2 == 1, H3K27me3_1 ==1))
nrow(filter(c42b_micro_3, NDR_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_3, NDR_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_3, NDR_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_3, NDR_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_3, H3K27me3_1 == 1, H3K27me3_2 ==1))
nrow(filter(c42b_micro_3, H3K27me3_1 == 1, H3K9me3_2 ==1)) + nrow(filter(c42b_micro_3, H3K27me3_2 == 1, H3K9me3_1 ==1))
nrow(filter(c42b_micro_3, H3K27me3_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_3, H3K27me3_2 == 1, None_1 ==1))


nrow(filter(c42b_micro_3, H3K9me3_1 == 1, H3K9me3_2 ==1))
nrow(filter(c42b_micro_3, H3K9me3_1 == 1, None_2 ==1)) + nrow(filter(c42b_micro_3, H3K9me3_2 == 1, None_1 ==1))

nrow(filter(c42b_micro_3, None_1 == 1, None_2 ==1))

c42b_micro_3_pp <- c42b_micro_3 %>% filter(Promoter_1 == 1, Promoter_2 ==1)
c42b_micro_3_pe1 <- c42b_micro_3 %>% filter(Promoter_1 == 1, Enhancer_2 ==1)
c42b_micro_3_pe2 <- c42b_micro_3 %>% filter(Promoter_2 == 1, Enhancer_1 ==1)
c42b_micro_3_pe <- rbind(c42b_micro_3_pe1,c42b_micro_3_pe2)

c42b_micro_3_pi1 <- c42b_micro_3 %>% filter(Promoter_1 == 1, Insulator_2 ==1)
c42b_micro_3_pi2 <- c42b_micro_3 %>% filter(Promoter_2 == 1, Insulator_1 ==1)
c42b_micro_3_pi <- rbind(c42b_micro_3_pi1,c42b_micro_3_pi2)

c42b_micro_3_pn1 <- c42b_micro_3 %>% filter(Promoter_1 == 1, NDR_2 ==1)
c42b_micro_3_pn2 <- c42b_micro_3 %>% filter(Promoter_2 == 1, NDR_1 ==1)
c42b_micro_3_pn <- rbind(c42b_micro_3_pn1,c42b_micro_3_pn2)

c42b_micro_3_ph271 <- c42b_micro_3 %>% filter(Promoter_1 == 1, H3K27me3_2 ==1)
c42b_micro_3_ph272 <- c42b_micro_3 %>% filter(Promoter_2 == 1, H3K27me3_1 ==1)
c42b_micro_3_ph27 <- rbind(c42b_micro_3_ph271,c42b_micro_3_ph272)

c42b_micro_3_ph91 <- c42b_micro_3 %>% filter(Promoter_1 == 1, H3K9me3_2 ==1)
c42b_micro_3_ph92 <- c42b_micro_3 %>% filter(Promoter_2 == 1, H3K9me3_1 ==1)
c42b_micro_3_ph9 <- rbind(c42b_micro_3_ph91,c42b_micro_3_ph92)

c42b_micro_3_px1 <- c42b_micro_3 %>% filter(Promoter_1 == 1, None_2 ==1)
c42b_micro_3_px2 <- c42b_micro_3 %>% filter(Promoter_2 == 1, None_1 ==1)
c42b_micro_3_px <- rbind(c42b_micro_3_px1,c42b_micro_3_px2)

c42b_micro_3_ee <- c42b_micro_3 %>% filter(Enhancer_1 == 1, Enhancer_2 ==1)

c42b_micro_3_ei1 <- c42b_micro_3 %>% filter(Enhancer_1 == 1, Insulator_2 ==1)
c42b_micro_3_ei2 <- c42b_micro_3 %>% filter(Enhancer_2 == 1, Insulator_1 ==1)
c42b_micro_3_ei <- rbind(c42b_micro_3_ei1,c42b_micro_3_ei2)

c42b_micro_3_en1 <- c42b_micro_3 %>% filter(Enhancer_1 == 1, NDR_2 ==1)
c42b_micro_3_en2 <- c42b_micro_3 %>% filter(Enhancer_2 == 1, NDR_1 ==1)
c42b_micro_3_en <- rbind(c42b_micro_3_en1,c42b_micro_3_en2)

c42b_micro_3_eh271 <- c42b_micro_3 %>% filter(Enhancer_1 == 1, H3K27me3_2 ==1)
c42b_micro_3_eh272 <- c42b_micro_3 %>% filter(Enhancer_2 == 1, H3K27me3_1 ==1)
c42b_micro_3_eh27 <- rbind(c42b_micro_3_eh271,c42b_micro_3_eh272)

c42b_micro_3_eh91 <- c42b_micro_3 %>% filter(Enhancer_1 == 1, H3K9me3_2 ==1)
c42b_micro_3_eh92 <- c42b_micro_3 %>% filter(Enhancer_2 == 1, H3K9me3_1 ==1)
c42b_micro_3_eh9 <- rbind(c42b_micro_3_eh91,c42b_micro_3_eh92)

c42b_micro_3_ex1 <- c42b_micro_3 %>% filter(Enhancer_1 == 1, None_2 ==1)
c42b_micro_3_ex2 <- c42b_micro_3 %>% filter(Enhancer_2 == 1, None_1 ==1)
c42b_micro_3_ex <- rbind(c42b_micro_3_ex1,c42b_micro_3_ex2)

c42b_micro_3_ii <- c42b_micro_3 %>% filter(Insulator_1 == 1, Insulator_2 ==1)

c42b_micro_3_in1 <- c42b_micro_3 %>% filter(Insulator_1 == 1, NDR_2 ==1)
c42b_micro_3_in2 <- c42b_micro_3 %>% filter(Insulator_2 == 1, NDR_1 ==1)
c42b_micro_3_in <- rbind(c42b_micro_3_in1,c42b_micro_3_in2)

c42b_micro_3_ih271 <- c42b_micro_3 %>% filter(Insulator_1 == 1, H3K27me3_2 ==1)
c42b_micro_3_ih272 <- c42b_micro_3 %>% filter(Insulator_2 == 1, H3K27me3_1 ==1)
c42b_micro_3_ih27 <- rbind(c42b_micro_3_ih271,c42b_micro_3_ih272)

c42b_micro_3_ih91 <- c42b_micro_3 %>% filter(Insulator_1 == 1, H3K9me3_2 ==1)
c42b_micro_3_ih92 <- c42b_micro_3 %>% filter(Insulator_2 == 1, H3K9me3_1 ==1)
c42b_micro_3_ih9 <- rbind(c42b_micro_3_ih91,c42b_micro_3_ih92)

c42b_micro_3_ix1 <- c42b_micro_3 %>% filter(Insulator_1 == 1, None_2 ==1)
c42b_micro_3_ix2 <- c42b_micro_3 %>% filter(Insulator_2 == 1, None_1 ==1)
c42b_micro_3_ix <- rbind(c42b_micro_3_ix1,c42b_micro_3_ix2)

c42b_micro_3_nn <- c42b_micro_3 %>% filter(NDR_1 == 1, NDR_2 ==1)


c42b_micro_3_nh271 <- c42b_micro_3 %>% filter(NDR_1 == 1, H3K27me3_2 ==1)
c42b_micro_3_nh272 <- c42b_micro_3 %>% filter(NDR_2 == 1, H3K27me3_1 ==1)
c42b_micro_3_nh27 <- rbind(c42b_micro_3_nh271,c42b_micro_3_nh272)

c42b_micro_3_nh91 <- c42b_micro_3 %>% filter(NDR_1 == 1, H3K9me3_2 ==1)
c42b_micro_3_nh92 <- c42b_micro_3 %>% filter(NDR_2 == 1, H3K9me3_1 ==1)
c42b_micro_3_nh9 <- rbind(c42b_micro_3_nh91,c42b_micro_3_nh92)

c42b_micro_3_nx1 <- c42b_micro_3 %>% filter(NDR_1 == 1, None_2 ==1)
c42b_micro_3_nx2 <- c42b_micro_3 %>% filter(NDR_2 == 1, None_1 ==1)
c42b_micro_3_nx <- rbind(c42b_micro_3_nx1,c42b_micro_3_nx2)


c42b_micro_3_h27h27 <- c42b_micro_3 %>% filter(H3K27me3_1 == 1, H3K27me3_2 ==1)

c42b_micro_3_h27h91 <- c42b_micro_3 %>% filter(H3K27me3_1 == 1, H3K9me3_2 ==1)
c42b_micro_3_h27h92 <- c42b_micro_3 %>% filter(H3K27me3_2 == 1, H3K9me3_1 ==1)
c42b_micro_3_h27h9 <- rbind(c42b_micro_3_h27h91,c42b_micro_3_h27h92)

c42b_micro_3_h27x1 <- c42b_micro_3 %>% filter(H3K27me3_1 == 1, None_2 ==1)
c42b_micro_3_h27x2 <- c42b_micro_3 %>% filter(H3K27me3_2 == 1, None_1 ==1)
c42b_micro_3_h27x <- rbind(c42b_micro_3_h27x1,c42b_micro_3_h27x2)

c42b_micro_3_h9h9 <- c42b_micro_3 %>% filter(H3K9me3_1 == 1, H3K9me3_2 ==1)

c42b_micro_3_h9x1 <- c42b_micro_3 %>% filter(H3K9me3_1 == 1, None_2 ==1)
c42b_micro_3_h9x2 <- c42b_micro_3 %>% filter(H3K9me3_2 == 1, None_1 ==1)
c42b_micro_3_h9x <- rbind(c42b_micro_3_h9x1,c42b_micro_3_h9x2)

c42b_micro_3_xx <- c42b_micro_3 %>% filter(None_1 == 1, None_2 ==1)

write.table(c42b_micro_3_pp, "Micro_C_3_Billion_5kb_Promoter_Promoter_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_pe, "Micro_C_3_Billion_5kb_Promoter_Enhancer_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_pi, "Micro_C_3_Billion_5kb_Promoter_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_pn, "Micro_C_3_Billion_5kb_Promoter_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_ph27, "Micro_C_3_Billion_5kb_Promoter_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_ph9, "Micro_C_3_Billion_5kb_Promoter_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_px, "Micro_C_3_Billion_5kb_Promoter_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_3_ee, "Micro_C_3_Billion_5kb_Enhancer_Enhancer_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_ei, "Micro_C_3_Billion_5kb_Enhancer_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_en, "Micro_C_3_Billion_5kb_Enhancer_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_eh27, "Micro_C_3_Billion_5kb_Enhancer_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_eh9, "Micro_C_3_Billion_5kb_Enhancer_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_ex, "Micro_C_3_Billion_5kb_Enhancer_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_3_ii, "Micro_C_3_Billion_5kb_Insulator_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_in, "Micro_C_3_Billion_5kb_Insulator_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_ih27, "Micro_C_3_Billion_5kb_Insulator_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_ih9, "Micro_C_3_Billion_5kb_Insulator_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_ix, "Micro_C_3_Billion_5kb_Insulator_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_3_nn, "Micro_C_3_Billion_5kb_NDR_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_nh27, "Micro_C_3_Billion_5kb_NDR_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_nh9, "Micro_C_3_Billion_5kb_NDR_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_nx, "Micro_C_3_Billion_5kb_NDR_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_3_h27h27, "Micro_C_3_Billion_5kb_H3K27me3_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_h27h9, "Micro_C_3_Billion_5kb_H3K27me3_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_h27x, "Micro_C_3_Billion_5kb_H3K27me3_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_3_h9h9, "Micro_C_3_Billion_5kb_H3K9me3_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_3_h9x, "Micro_C_3_Billion_5kb_H3K9me3_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_3_xx, "Micro_C_3_Billion_5kb_None_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

c42b_micro_2_pp <- c42b_micro_2 %>% filter(Promoter_1 == 1, Promoter_2 ==1)
c42b_micro_2_pe1 <- c42b_micro_2 %>% filter(Promoter_1 == 1, Enhancer_2 ==1)
c42b_micro_2_pe2 <- c42b_micro_2 %>% filter(Promoter_2 == 1, Enhancer_1 ==1)
c42b_micro_2_pe <- rbind(c42b_micro_2_pe1,c42b_micro_2_pe2)

c42b_micro_2_pi1 <- c42b_micro_2 %>% filter(Promoter_1 == 1, Insulator_2 ==1)
c42b_micro_2_pi2 <- c42b_micro_2 %>% filter(Promoter_2 == 1, Insulator_1 ==1)
c42b_micro_2_pi <- rbind(c42b_micro_2_pi1,c42b_micro_2_pi2)

c42b_micro_2_pn1 <- c42b_micro_2 %>% filter(Promoter_1 == 1, NDR_2 ==1)
c42b_micro_2_pn2 <- c42b_micro_2 %>% filter(Promoter_2 == 1, NDR_1 ==1)
c42b_micro_2_pn <- rbind(c42b_micro_2_pn1,c42b_micro_2_pn2)

c42b_micro_2_ph271 <- c42b_micro_2 %>% filter(Promoter_1 == 1, H3K27me3_2 ==1)
c42b_micro_2_ph272 <- c42b_micro_2 %>% filter(Promoter_2 == 1, H3K27me3_1 ==1)
c42b_micro_2_ph27 <- rbind(c42b_micro_2_ph271,c42b_micro_2_ph272)

c42b_micro_2_ph91 <- c42b_micro_2 %>% filter(Promoter_1 == 1, H3K9me3_2 ==1)
c42b_micro_2_ph92 <- c42b_micro_2 %>% filter(Promoter_2 == 1, H3K9me3_1 ==1)
c42b_micro_2_ph9 <- rbind(c42b_micro_2_ph91,c42b_micro_2_ph92)

c42b_micro_2_px1 <- c42b_micro_2 %>% filter(Promoter_1 == 1, None_2 ==1)
c42b_micro_2_px2 <- c42b_micro_2 %>% filter(Promoter_2 == 1, None_1 ==1)
c42b_micro_2_px <- rbind(c42b_micro_2_px1,c42b_micro_2_px2)

c42b_micro_2_ee <- c42b_micro_2 %>% filter(Enhancer_1 == 1, Enhancer_2 ==1)

c42b_micro_2_ei1 <- c42b_micro_2 %>% filter(Enhancer_1 == 1, Insulator_2 ==1)
c42b_micro_2_ei2 <- c42b_micro_2 %>% filter(Enhancer_2 == 1, Insulator_1 ==1)
c42b_micro_2_ei <- rbind(c42b_micro_2_ei1,c42b_micro_2_ei2)

c42b_micro_2_en1 <- c42b_micro_2 %>% filter(Enhancer_1 == 1, NDR_2 ==1)
c42b_micro_2_en2 <- c42b_micro_2 %>% filter(Enhancer_2 == 1, NDR_1 ==1)
c42b_micro_2_en <- rbind(c42b_micro_2_en1,c42b_micro_2_en2)

c42b_micro_2_eh271 <- c42b_micro_2 %>% filter(Enhancer_1 == 1, H3K27me3_2 ==1)
c42b_micro_2_eh272 <- c42b_micro_2 %>% filter(Enhancer_2 == 1, H3K27me3_1 ==1)
c42b_micro_2_eh27 <- rbind(c42b_micro_2_eh271,c42b_micro_2_eh272)

c42b_micro_2_eh91 <- c42b_micro_2 %>% filter(Enhancer_1 == 1, H3K9me3_2 ==1)
c42b_micro_2_eh92 <- c42b_micro_2 %>% filter(Enhancer_2 == 1, H3K9me3_1 ==1)
c42b_micro_2_eh9 <- rbind(c42b_micro_2_eh91,c42b_micro_2_eh92)

c42b_micro_2_ex1 <- c42b_micro_2 %>% filter(Enhancer_1 == 1, None_2 ==1)
c42b_micro_2_ex2 <- c42b_micro_2 %>% filter(Enhancer_2 == 1, None_1 ==1)
c42b_micro_2_ex <- rbind(c42b_micro_2_ex1,c42b_micro_2_ex2)

c42b_micro_2_ii <- c42b_micro_2 %>% filter(Insulator_1 == 1, Insulator_2 ==1)

c42b_micro_2_in1 <- c42b_micro_2 %>% filter(Insulator_1 == 1, NDR_2 ==1)
c42b_micro_2_in2 <- c42b_micro_2 %>% filter(Insulator_2 == 1, NDR_1 ==1)
c42b_micro_2_in <- rbind(c42b_micro_2_in1,c42b_micro_2_in2)

c42b_micro_2_ih271 <- c42b_micro_2 %>% filter(Insulator_1 == 1, H3K27me3_2 ==1)
c42b_micro_2_ih272 <- c42b_micro_2 %>% filter(Insulator_2 == 1, H3K27me3_1 ==1)
c42b_micro_2_ih27 <- rbind(c42b_micro_2_ih271,c42b_micro_2_ih272)

c42b_micro_2_ih91 <- c42b_micro_2 %>% filter(Insulator_1 == 1, H3K9me3_2 ==1)
c42b_micro_2_ih92 <- c42b_micro_2 %>% filter(Insulator_2 == 1, H3K9me3_1 ==1)
c42b_micro_2_ih9 <- rbind(c42b_micro_2_ih91,c42b_micro_2_ih92)

c42b_micro_2_ix1 <- c42b_micro_2 %>% filter(Insulator_1 == 1, None_2 ==1)
c42b_micro_2_ix2 <- c42b_micro_2 %>% filter(Insulator_2 == 1, None_1 ==1)
c42b_micro_2_ix <- rbind(c42b_micro_2_ix1,c42b_micro_2_ix2)

c42b_micro_2_nn <- c42b_micro_2 %>% filter(NDR_1 == 1, NDR_2 ==1)


c42b_micro_2_nh271 <- c42b_micro_2 %>% filter(NDR_1 == 1, H3K27me3_2 ==1)
c42b_micro_2_nh272 <- c42b_micro_2 %>% filter(NDR_2 == 1, H3K27me3_1 ==1)
c42b_micro_2_nh27 <- rbind(c42b_micro_2_nh271,c42b_micro_2_nh272)

c42b_micro_2_nh91 <- c42b_micro_2 %>% filter(NDR_1 == 1, H3K9me3_2 ==1)
c42b_micro_2_nh92 <- c42b_micro_2 %>% filter(NDR_2 == 1, H3K9me3_1 ==1)
c42b_micro_2_nh9 <- rbind(c42b_micro_2_nh91,c42b_micro_2_nh92)

c42b_micro_2_nx1 <- c42b_micro_2 %>% filter(NDR_1 == 1, None_2 ==1)
c42b_micro_2_nx2 <- c42b_micro_2 %>% filter(NDR_2 == 1, None_1 ==1)
c42b_micro_2_nx <- rbind(c42b_micro_2_nx1,c42b_micro_2_nx2)


c42b_micro_2_h27h27 <- c42b_micro_2 %>% filter(H3K27me3_1 == 1, H3K27me3_2 ==1)

c42b_micro_2_h27h91 <- c42b_micro_2 %>% filter(H3K27me3_1 == 1, H3K9me3_2 ==1)
c42b_micro_2_h27h92 <- c42b_micro_2 %>% filter(H3K27me3_2 == 1, H3K9me3_1 ==1)
c42b_micro_2_h27h9 <- rbind(c42b_micro_2_h27h91,c42b_micro_2_h27h92)

c42b_micro_2_h27x1 <- c42b_micro_2 %>% filter(H3K27me3_1 == 1, None_2 ==1)
c42b_micro_2_h27x2 <- c42b_micro_2 %>% filter(H3K27me3_2 == 1, None_1 ==1)
c42b_micro_2_h27x <- rbind(c42b_micro_2_h27x1,c42b_micro_2_h27x2)

c42b_micro_2_h9h9 <- c42b_micro_2 %>% filter(H3K9me3_1 == 1, H3K9me3_2 ==1)

c42b_micro_2_h9x1 <- c42b_micro_2 %>% filter(H3K9me3_1 == 1, None_2 ==1)
c42b_micro_2_h9x2 <- c42b_micro_2 %>% filter(H3K9me3_2 == 1, None_1 ==1)
c42b_micro_2_h9x <- rbind(c42b_micro_2_h9x1,c42b_micro_2_h9x2)

c42b_micro_2_xx <- c42b_micro_2 %>% filter(None_1 == 1, None_2 ==1)

write.table(c42b_micro_2_pp, "Micro_C_2_Billion_5kb_Promoter_Promoter_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_pe, "Micro_C_2_Billion_5kb_Promoter_Enhancer_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_pi, "Micro_C_2_Billion_5kb_Promoter_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_pn, "Micro_C_2_Billion_5kb_Promoter_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_ph27, "Micro_C_2_Billion_5kb_Promoter_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_ph9, "Micro_C_2_Billion_5kb_Promoter_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_px, "Micro_C_2_Billion_5kb_Promoter_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_2_ee, "Micro_C_2_Billion_5kb_Enhancer_Enhancer_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_ei, "Micro_C_2_Billion_5kb_Enhancer_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_en, "Micro_C_2_Billion_5kb_Enhancer_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_eh27, "Micro_C_2_Billion_5kb_Enhancer_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_eh9, "Micro_C_2_Billion_5kb_Enhancer_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_ex, "Micro_C_2_Billion_5kb_Enhancer_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_2_ii, "Micro_C_2_Billion_5kb_Insulator_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_in, "Micro_C_2_Billion_5kb_Insulator_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_ih27, "Micro_C_2_Billion_5kb_Insulator_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_ih9, "Micro_C_2_Billion_5kb_Insulator_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_ix, "Micro_C_2_Billion_5kb_Insulator_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_2_nn, "Micro_C_2_Billion_5kb_NDR_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_nh27, "Micro_C_2_Billion_5kb_NDR_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_nh9, "Micro_C_2_Billion_5kb_NDR_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_nx, "Micro_C_2_Billion_5kb_NDR_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_2_h27h27, "Micro_C_2_Billion_5kb_H3K27me3_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_h27h9, "Micro_C_2_Billion_5kb_H3K27me3_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_h27x, "Micro_C_2_Billion_5kb_H3K27me3_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_2_h9h9, "Micro_C_2_Billion_5kb_H3K9me3_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_2_h9x, "Micro_C_2_Billion_5kb_H3K9me3_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_2_xx, "Micro_C_2_Billion_5kb_None_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

c42b_micro_1_pp <- c42b_micro_1 %>% filter(Promoter_1 == 1, Promoter_2 ==1)
c42b_micro_1_pe1 <- c42b_micro_1 %>% filter(Promoter_1 == 1, Enhancer_2 ==1)
c42b_micro_1_pe2 <- c42b_micro_1 %>% filter(Promoter_2 == 1, Enhancer_1 ==1)
c42b_micro_1_pe <- rbind(c42b_micro_1_pe1,c42b_micro_1_pe2)

c42b_micro_1_pi1 <- c42b_micro_1 %>% filter(Promoter_1 == 1, Insulator_2 ==1)
c42b_micro_1_pi2 <- c42b_micro_1 %>% filter(Promoter_2 == 1, Insulator_1 ==1)
c42b_micro_1_pi <- rbind(c42b_micro_1_pi1,c42b_micro_1_pi2)

c42b_micro_1_pn1 <- c42b_micro_1 %>% filter(Promoter_1 == 1, NDR_2 ==1)
c42b_micro_1_pn2 <- c42b_micro_1 %>% filter(Promoter_2 == 1, NDR_1 ==1)
c42b_micro_1_pn <- rbind(c42b_micro_1_pn1,c42b_micro_1_pn2)

c42b_micro_1_ph271 <- c42b_micro_1 %>% filter(Promoter_1 == 1, H3K27me3_2 ==1)
c42b_micro_1_ph272 <- c42b_micro_1 %>% filter(Promoter_2 == 1, H3K27me3_1 ==1)
c42b_micro_1_ph27 <- rbind(c42b_micro_1_ph271,c42b_micro_1_ph272)

c42b_micro_1_ph91 <- c42b_micro_1 %>% filter(Promoter_1 == 1, H3K9me3_2 ==1)
c42b_micro_1_ph92 <- c42b_micro_1 %>% filter(Promoter_2 == 1, H3K9me3_1 ==1)
c42b_micro_1_ph9 <- rbind(c42b_micro_1_ph91,c42b_micro_1_ph92)

c42b_micro_1_px1 <- c42b_micro_1 %>% filter(Promoter_1 == 1, None_2 ==1)
c42b_micro_1_px2 <- c42b_micro_1 %>% filter(Promoter_2 == 1, None_1 ==1)
c42b_micro_1_px <- rbind(c42b_micro_1_px1,c42b_micro_1_px2)

c42b_micro_1_ee <- c42b_micro_1 %>% filter(Enhancer_1 == 1, Enhancer_2 ==1)

c42b_micro_1_ei1 <- c42b_micro_1 %>% filter(Enhancer_1 == 1, Insulator_2 ==1)
c42b_micro_1_ei2 <- c42b_micro_1 %>% filter(Enhancer_2 == 1, Insulator_1 ==1)
c42b_micro_1_ei <- rbind(c42b_micro_1_ei1,c42b_micro_1_ei2)

c42b_micro_1_en1 <- c42b_micro_1 %>% filter(Enhancer_1 == 1, NDR_2 ==1)
c42b_micro_1_en2 <- c42b_micro_1 %>% filter(Enhancer_2 == 1, NDR_1 ==1)
c42b_micro_1_en <- rbind(c42b_micro_1_en1,c42b_micro_1_en2)

c42b_micro_1_eh271 <- c42b_micro_1 %>% filter(Enhancer_1 == 1, H3K27me3_2 ==1)
c42b_micro_1_eh272 <- c42b_micro_1 %>% filter(Enhancer_2 == 1, H3K27me3_1 ==1)
c42b_micro_1_eh27 <- rbind(c42b_micro_1_eh271,c42b_micro_1_eh272)

c42b_micro_1_eh91 <- c42b_micro_1 %>% filter(Enhancer_1 == 1, H3K9me3_2 ==1)
c42b_micro_1_eh92 <- c42b_micro_1 %>% filter(Enhancer_2 == 1, H3K9me3_1 ==1)
c42b_micro_1_eh9 <- rbind(c42b_micro_1_eh91,c42b_micro_1_eh92)

c42b_micro_1_ex1 <- c42b_micro_1 %>% filter(Enhancer_1 == 1, None_2 ==1)
c42b_micro_1_ex2 <- c42b_micro_1 %>% filter(Enhancer_2 == 1, None_1 ==1)
c42b_micro_1_ex <- rbind(c42b_micro_1_ex1,c42b_micro_1_ex2)

c42b_micro_1_ii <- c42b_micro_1 %>% filter(Insulator_1 == 1, Insulator_2 ==1)

c42b_micro_1_in1 <- c42b_micro_1 %>% filter(Insulator_1 == 1, NDR_2 ==1)
c42b_micro_1_in2 <- c42b_micro_1 %>% filter(Insulator_2 == 1, NDR_1 ==1)
c42b_micro_1_in <- rbind(c42b_micro_1_in1,c42b_micro_1_in2)

c42b_micro_1_ih271 <- c42b_micro_1 %>% filter(Insulator_1 == 1, H3K27me3_2 ==1)
c42b_micro_1_ih272 <- c42b_micro_1 %>% filter(Insulator_2 == 1, H3K27me3_1 ==1)
c42b_micro_1_ih27 <- rbind(c42b_micro_1_ih271,c42b_micro_1_ih272)

c42b_micro_1_ih91 <- c42b_micro_1 %>% filter(Insulator_1 == 1, H3K9me3_2 ==1)
c42b_micro_1_ih92 <- c42b_micro_1 %>% filter(Insulator_2 == 1, H3K9me3_1 ==1)
c42b_micro_1_ih9 <- rbind(c42b_micro_1_ih91,c42b_micro_1_ih92)

c42b_micro_1_ix1 <- c42b_micro_1 %>% filter(Insulator_1 == 1, None_2 ==1)
c42b_micro_1_ix2 <- c42b_micro_1 %>% filter(Insulator_2 == 1, None_1 ==1)
c42b_micro_1_ix <- rbind(c42b_micro_1_ix1,c42b_micro_1_ix2)

c42b_micro_1_nn <- c42b_micro_1 %>% filter(NDR_1 == 1, NDR_2 ==1)


c42b_micro_1_nh271 <- c42b_micro_1 %>% filter(NDR_1 == 1, H3K27me3_2 ==1)
c42b_micro_1_nh272 <- c42b_micro_1 %>% filter(NDR_2 == 1, H3K27me3_1 ==1)
c42b_micro_1_nh27 <- rbind(c42b_micro_1_nh271,c42b_micro_1_nh272)

c42b_micro_1_nh91 <- c42b_micro_1 %>% filter(NDR_1 == 1, H3K9me3_2 ==1)
c42b_micro_1_nh92 <- c42b_micro_1 %>% filter(NDR_2 == 1, H3K9me3_1 ==1)
c42b_micro_1_nh9 <- rbind(c42b_micro_1_nh91,c42b_micro_1_nh92)

c42b_micro_1_nx1 <- c42b_micro_1 %>% filter(NDR_1 == 1, None_2 ==1)
c42b_micro_1_nx2 <- c42b_micro_1 %>% filter(NDR_2 == 1, None_1 ==1)
c42b_micro_1_nx <- rbind(c42b_micro_1_nx1,c42b_micro_1_nx2)


c42b_micro_1_h27h27 <- c42b_micro_1 %>% filter(H3K27me3_1 == 1, H3K27me3_2 ==1)

c42b_micro_1_h27h91 <- c42b_micro_1 %>% filter(H3K27me3_1 == 1, H3K9me3_2 ==1)
c42b_micro_1_h27h92 <- c42b_micro_1 %>% filter(H3K27me3_2 == 1, H3K9me3_1 ==1)
c42b_micro_1_h27h9 <- rbind(c42b_micro_1_h27h91,c42b_micro_1_h27h92)

c42b_micro_1_h27x1 <- c42b_micro_1 %>% filter(H3K27me3_1 == 1, None_2 ==1)
c42b_micro_1_h27x2 <- c42b_micro_1 %>% filter(H3K27me3_2 == 1, None_1 ==1)
c42b_micro_1_h27x <- rbind(c42b_micro_1_h27x1,c42b_micro_1_h27x2)

c42b_micro_1_h9h9 <- c42b_micro_1 %>% filter(H3K9me3_1 == 1, H3K9me3_2 ==1)

c42b_micro_1_h9x1 <- c42b_micro_1 %>% filter(H3K9me3_1 == 1, None_2 ==1)
c42b_micro_1_h9x2 <- c42b_micro_1 %>% filter(H3K9me3_2 == 1, None_1 ==1)
c42b_micro_1_h9x <- rbind(c42b_micro_1_h9x1,c42b_micro_1_h9x2)

c42b_micro_1_xx <- c42b_micro_1 %>% filter(None_1 == 1, None_2 ==1)

write.table(c42b_micro_1_pp, "Micro_C_1_Billion_5kb_Promoter_Promoter_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_pe, "Micro_C_1_Billion_5kb_Promoter_Enhancer_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_pi, "Micro_C_1_Billion_5kb_Promoter_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_pn, "Micro_C_1_Billion_5kb_Promoter_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_ph27, "Micro_C_1_Billion_5kb_Promoter_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_ph9, "Micro_C_1_Billion_5kb_Promoter_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_px, "Micro_C_1_Billion_5kb_Promoter_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_1_ee, "Micro_C_1_Billion_5kb_Enhancer_Enhancer_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_ei, "Micro_C_1_Billion_5kb_Enhancer_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_en, "Micro_C_1_Billion_5kb_Enhancer_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_eh27, "Micro_C_1_Billion_5kb_Enhancer_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_eh9, "Micro_C_1_Billion_5kb_Enhancer_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_ex, "Micro_C_1_Billion_5kb_Enhancer_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_1_ii, "Micro_C_1_Billion_5kb_Insulator_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_in, "Micro_C_1_Billion_5kb_Insulator_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_ih27, "Micro_C_1_Billion_5kb_Insulator_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_ih9, "Micro_C_1_Billion_5kb_Insulator_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_ix, "Micro_C_1_Billion_5kb_Insulator_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_1_nn, "Micro_C_1_Billion_5kb_NDR_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_nh27, "Micro_C_1_Billion_5kb_NDR_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_nh9, "Micro_C_1_Billion_5kb_NDR_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_nx, "Micro_C_1_Billion_5kb_NDR_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_1_h27h27, "Micro_C_1_Billion_5kb_H3K27me3_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_h27h9, "Micro_C_1_Billion_5kb_H3K27me3_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_h27x, "Micro_C_1_Billion_5kb_H3K27me3_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_1_h9h9, "Micro_C_1_Billion_5kb_H3K9me3_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_micro_1_h9x, "Micro_C_1_Billion_5kb_H3K9me3_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_micro_1_xx, "Micro_C_1_Billion_5kb_None_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

c42b_hic_pp <- c42b_hic %>% filter(Promoter_1 == 1, Promoter_2 ==1)
c42b_hic_pe1 <- c42b_hic %>% filter(Promoter_1 == 1, Enhancer_2 ==1)
c42b_hic_pe2 <- c42b_hic %>% filter(Promoter_2 == 1, Enhancer_1 ==1)
c42b_hic_pe <- rbind(c42b_hic_pe1,c42b_hic_pe2)

c42b_hic_pi1 <- c42b_hic %>% filter(Promoter_1 == 1, Insulator_2 ==1)
c42b_hic_pi2 <- c42b_hic %>% filter(Promoter_2 == 1, Insulator_1 ==1)
c42b_hic_pi <- rbind(c42b_hic_pi1,c42b_hic_pi2)

c42b_hic_pn1 <- c42b_hic %>% filter(Promoter_1 == 1, NDR_2 ==1)
c42b_hic_pn2 <- c42b_hic %>% filter(Promoter_2 == 1, NDR_1 ==1)
c42b_hic_pn <- rbind(c42b_hic_pn1,c42b_hic_pn2)

c42b_hic_ph271 <- c42b_hic %>% filter(Promoter_1 == 1, H3K27me3_2 ==1)
c42b_hic_ph272 <- c42b_hic %>% filter(Promoter_2 == 1, H3K27me3_1 ==1)
c42b_hic_ph27 <- rbind(c42b_hic_ph271,c42b_hic_ph272)

c42b_hic_ph91 <- c42b_hic %>% filter(Promoter_1 == 1, H3K9me3_2 ==1)
c42b_hic_ph92 <- c42b_hic %>% filter(Promoter_2 == 1, H3K9me3_1 ==1)
c42b_hic_ph9 <- rbind(c42b_hic_ph91,c42b_hic_ph92)

c42b_hic_px1 <- c42b_hic %>% filter(Promoter_1 == 1, None_2 ==1)
c42b_hic_px2 <- c42b_hic %>% filter(Promoter_2 == 1, None_1 ==1)
c42b_hic_px <- rbind(c42b_hic_px1,c42b_hic_px2)

c42b_hic_ee <- c42b_hic %>% filter(Enhancer_1 == 1, Enhancer_2 ==1)

c42b_hic_ei1 <- c42b_hic %>% filter(Enhancer_1 == 1, Insulator_2 ==1)
c42b_hic_ei2 <- c42b_hic %>% filter(Enhancer_2 == 1, Insulator_1 ==1)
c42b_hic_ei <- rbind(c42b_hic_ei1,c42b_hic_ei2)

c42b_hic_en1 <- c42b_hic %>% filter(Enhancer_1 == 1, NDR_2 ==1)
c42b_hic_en2 <- c42b_hic %>% filter(Enhancer_2 == 1, NDR_1 ==1)
c42b_hic_en <- rbind(c42b_hic_en1,c42b_hic_en2)

c42b_hic_eh271 <- c42b_hic %>% filter(Enhancer_1 == 1, H3K27me3_2 ==1)
c42b_hic_eh272 <- c42b_hic %>% filter(Enhancer_2 == 1, H3K27me3_1 ==1)
c42b_hic_eh27 <- rbind(c42b_hic_eh271,c42b_hic_eh272)

c42b_hic_eh91 <- c42b_hic %>% filter(Enhancer_1 == 1, H3K9me3_2 ==1)
c42b_hic_eh92 <- c42b_hic %>% filter(Enhancer_2 == 1, H3K9me3_1 ==1)
c42b_hic_eh9 <- rbind(c42b_hic_eh91,c42b_hic_eh92)

c42b_hic_ex1 <- c42b_hic %>% filter(Enhancer_1 == 1, None_2 ==1)
c42b_hic_ex2 <- c42b_hic %>% filter(Enhancer_2 == 1, None_1 ==1)
c42b_hic_ex <- rbind(c42b_hic_ex1,c42b_hic_ex2)

c42b_hic_ii <- c42b_hic %>% filter(Insulator_1 == 1, Insulator_2 ==1)

c42b_hic_in1 <- c42b_hic %>% filter(Insulator_1 == 1, NDR_2 ==1)
c42b_hic_in2 <- c42b_hic %>% filter(Insulator_2 == 1, NDR_1 ==1)
c42b_hic_in <- rbind(c42b_hic_in1,c42b_hic_in2)

c42b_hic_ih271 <- c42b_hic %>% filter(Insulator_1 == 1, H3K27me3_2 ==1)
c42b_hic_ih272 <- c42b_hic %>% filter(Insulator_2 == 1, H3K27me3_1 ==1)
c42b_hic_ih27 <- rbind(c42b_hic_ih271,c42b_hic_ih272)

c42b_hic_ih91 <- c42b_hic %>% filter(Insulator_1 == 1, H3K9me3_2 ==1)
c42b_hic_ih92 <- c42b_hic %>% filter(Insulator_2 == 1, H3K9me3_1 ==1)
c42b_hic_ih9 <- rbind(c42b_hic_ih91,c42b_hic_ih92)

c42b_hic_ix1 <- c42b_hic %>% filter(Insulator_1 == 1, None_2 ==1)
c42b_hic_ix2 <- c42b_hic %>% filter(Insulator_2 == 1, None_1 ==1)
c42b_hic_ix <- rbind(c42b_hic_ix1,c42b_hic_ix2)

c42b_hic_nn <- c42b_hic %>% filter(NDR_1 == 1, NDR_2 ==1)


c42b_hic_nh271 <- c42b_hic %>% filter(NDR_1 == 1, H3K27me3_2 ==1)
c42b_hic_nh272 <- c42b_hic %>% filter(NDR_2 == 1, H3K27me3_1 ==1)
c42b_hic_nh27 <- rbind(c42b_hic_nh271,c42b_hic_nh272)

c42b_hic_nh91 <- c42b_hic %>% filter(NDR_1 == 1, H3K9me3_2 ==1)
c42b_hic_nh92 <- c42b_hic %>% filter(NDR_2 == 1, H3K9me3_1 ==1)
c42b_hic_nh9 <- rbind(c42b_hic_nh91,c42b_hic_nh92)

c42b_hic_nx1 <- c42b_hic %>% filter(NDR_1 == 1, None_2 ==1)
c42b_hic_nx2 <- c42b_hic %>% filter(NDR_2 == 1, None_1 ==1)
c42b_hic_nx <- rbind(c42b_hic_nx1,c42b_hic_nx2)


c42b_hic_h27h27 <- c42b_hic %>% filter(H3K27me3_1 == 1, H3K27me3_2 ==1)

c42b_hic_h27h91 <- c42b_hic %>% filter(H3K27me3_1 == 1, H3K9me3_2 ==1)
c42b_hic_h27h92 <- c42b_hic %>% filter(H3K27me3_2 == 1, H3K9me3_1 ==1)
c42b_hic_h27h9 <- rbind(c42b_hic_h27h91,c42b_hic_h27h92)

c42b_hic_h27x1 <- c42b_hic %>% filter(H3K27me3_1 == 1, None_2 ==1)
c42b_hic_h27x2 <- c42b_hic %>% filter(H3K27me3_2 == 1, None_1 ==1)
c42b_hic_h27x <- rbind(c42b_hic_h27x1,c42b_hic_h27x2)

c42b_hic_h9h9 <- c42b_hic %>% filter(H3K9me3_1 == 1, H3K9me3_2 ==1)

c42b_hic_h9x1 <- c42b_hic %>% filter(H3K9me3_1 == 1, None_2 ==1)
c42b_hic_h9x2 <- c42b_hic %>% filter(H3K9me3_2 == 1, None_1 ==1)
c42b_hic_h9x <- rbind(c42b_hic_h9x1,c42b_hic_h9x2)

c42b_hic_xx <- c42b_hic %>% filter(None_1 == 1, None_2 ==1)

write.table(c42b_hic_pp, "HiC_Billion_5kb_Promoter_Promoter_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_pe, "HiC_Billion_5kb_Promoter_Enhancer_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_pi, "HiC_Billion_5kb_Promoter_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_pn, "HiC_Billion_5kb_Promoter_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_ph27, "HiC_Billion_5kb_Promoter_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_ph9, "HiC_Billion_5kb_Promoter_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_px, "HiC_Billion_5kb_Promoter_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_hic_ee, "HiC_Billion_5kb_Enhancer_Enhancer_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_ei, "HiC_Billion_5kb_Enhancer_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_en, "HiC_Billion_5kb_Enhancer_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_eh27, "HiC_Billion_5kb_Enhancer_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_eh9, "HiC_Billion_5kb_Enhancer_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_ex, "HiC_Billion_5kb_Enhancer_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_hic_ii, "HiC_Billion_5kb_Insulator_Insulator_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_in, "HiC_Billion_5kb_Insulator_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_ih27, "HiC_Billion_5kb_Insulator_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_ih9, "HiC_Billion_5kb_Insulator_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_ix, "HiC_Billion_5kb_Insulator_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_hic_nn, "HiC_Billion_5kb_NDR_NDR_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_nh27, "HiC_Billion_5kb_NDR_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_nh9, "HiC_Billion_5kb_NDR_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_nx, "HiC_Billion_5kb_NDR_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_hic_h27h27, "HiC_Billion_5kb_H3K27me3_H3K27me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_h27h9, "HiC_Billion_5kb_H3K27me3_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_h27x, "HiC_Billion_5kb_H3K27me3_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_hic_h9h9, "HiC_Billion_5kb_H3K9me3_H3K9me3_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(c42b_hic_h9x, "HiC_Billion_5kb_H3K9me3_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(c42b_hic_xx, "HiC_Billion_5kb_None_None_Loop.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")


