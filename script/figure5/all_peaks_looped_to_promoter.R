library(dplyr)
library(plyr)
library(fuzzyjoin)
library(ggplot2)
library(ggpubr)

c42b_micro_2 <- read.table("2_Billion_Mustache-5kb-all-loop.bedpe")
c42b_chicago <-  read.table("PC-UNI1945.5kb.bedpe")
c42b_hic <-  read.table("HiC_Mustache-5kb-all-loop.bedpe")
c42b_micro_1 <- read.table("1_Billion_Mustache-5kb-all-loop.bedpe")
c42b_micro_3 <- read.table("3_billion_low_Mustache-5kb-all-loop.bedpe")

c42b_micro_2$loopname <- paste("c42b_2_billion",rownames(c42b_micro_2),sep = "_")
c42b_micro_1$loopname <- paste("c42b_1_billion",rownames(c42b_micro_1),sep = "_")
c42b_chicago$loopname <- paste("c42b_chicago",rownames(c42b_chicago),sep = "_")
c42b_hic$loopname <- paste("c42b_hic",rownames(c42b_hic),sep = "_")
c42b_micro_3$loopname <- paste("c42b_3_billion",rownames(c42b_micro_3),sep = "_")

c42b_micro_3$V1 <- gsub("chr", "", c42b_micro_3$V1)
c42b_micro_3$V4 <- gsub("chr", "", c42b_micro_3$V4)
c42b_chicago$V1 <- gsub("chr", "", c42b_chicago$V1)
c42b_chicago$V4 <- gsub("chr", "", c42b_chicago$V4)

h3k27ac <- read.table("C42B_H3K27ac_no_tss.bed")
h3k4me3 <- read.table("C42B_TSS_FPKM_greater_than_0.5_with_enst.tsv")
ctcf <- read.table("C42B_CTCF_no_tss_no_enhancer.bed")
ndr <- read.table("C42B_NDR_without_features.bed")

ctcf$V1 <- gsub("chr","",ctcf$V1) 
h3k27ac$V1 <- gsub("chr","",h3k27ac$V1)
h3k4me3$V1 <- gsub("chr","",h3k4me3$V1)
ndr$V1 <- gsub("chr","",ndr$V1)

h3k4me3$peak_name <- paste("Peak",rownames(h3k4me3),sep ="_")

x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")


h3k27me3 <- read.table("H3K27me3_not_in_enhancer_insulator_promoter_ndr_region.bed")
h3k9me3 <- read.table("H3K9me3_not_in_enhancer_insulator_promoter_ndr_h3k27me3_region.bed")

ctcf$V1 <- gsub("chr","",ctcf$V1) 
h3k27ac$V1 <- gsub("chr","",h3k27ac$V1)
h3k4me3$V1 <- gsub("chr","",h3k4me3$V1)
ndr$V1 <- gsub("chr","",ndr$V1)
h3k27me3$V1 <- gsub("chr","",h3k27me3$V1)
h3k9me3$V1 <- gsub("chr","",h3k9me3$V1)

c42b_chicago$V1 <- gsub("chr","",c42b_chicago$V1)
c42b_chicago$V4 <- gsub("chr","",c42b_chicago$V4)

c42b_micro_3 <- c42b_micro_3 %>% mutate(Enhancer_1 = 0, Enhancer_2 = 0, Promoter_1 = 0, Promoter_2 = 0, Insulator_1 = 0, Insulator_2 = 0, NDR_1 = 0, NDR_2 = 0, H3K27me3_1 = 0, H3K27me3_2 = 0, H3K9me3_1 = 0, H3K9me3_2 = 0, None_1 = 1, None_2 = 1)




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



c42b_micro3_p1e2 <- c42b_micro_3 %>% filter(Promoter_1 == 1 & Enhancer_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_p2e1 <- c42b_micro_3 %>% filter(Promoter_2 == 1 & Enhancer_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_promoter1_with_enhancer <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p1e2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter1_with_enhancer <-rbind(c42b_micro_3_promoter1_with_enhancer,some_count)
}

c42b_micro3_p2e1 <- c42b_micro3_p2e1 %>% select(V4,V5,V6)
colnames(c42b_micro3_p2e1)[1] <- "V1"
colnames(c42b_micro3_p2e1)[2] <- "V2"
colnames(c42b_micro3_p2e1)[3] <- "V3"

c42b_micro_3_promoter2_with_enhancer <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p2e1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter2_with_enhancer <-rbind(c42b_micro_3_promoter2_with_enhancer,some_count)
}

c42b_pe <- rbind(c42b_micro_3_promoter1_with_enhancer, c42b_micro_3_promoter2_with_enhancer)
c42b_pe <- unique(c42b_pe)

write.table(c42b_pe, "Promoter_in_loop_with_enhancer.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")




c42b_micro3_p1i2 <- c42b_micro_3 %>% filter(Promoter_1 == 1 & Insulator_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_p2i1 <- c42b_micro_3 %>% filter(Promoter_2 == 1 & Insulator_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_promoter1_with_insulator <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p1i2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter1_with_insulator <-rbind(c42b_micro_3_promoter1_with_insulator,some_count)
}

c42b_micro3_p2i1 <- c42b_micro3_p2i1 %>% select(V4,V5,V6)
colnames(c42b_micro3_p2i1)[1] <- "V1"
colnames(c42b_micro3_p2i1)[2] <- "V2"
colnames(c42b_micro3_p2i1)[3] <- "V3"

c42b_micro_3_promoter2_with_insulator <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p2i1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter2_with_insulator <-rbind(c42b_micro_3_promoter2_with_insulator,some_count)
}

c42b_pi <- rbind(c42b_micro_3_promoter1_with_insulator, c42b_micro_3_promoter2_with_insulator)
c42b_pi <- unique(c42b_pi)

write.table(c42b_pi, "Promoter_in_loop_with_insulator.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")




c42b_micro3_p1e2 <- c42b_micro_3 %>% filter(Promoter_1 == 1 & NDR_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_p2e1 <- c42b_micro_3 %>% filter(Promoter_2 == 1 & NDR_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_promoter1_with_ndr <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p1e2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter1_with_ndr <-rbind(c42b_micro_3_promoter1_with_ndr,some_count)
}

c42b_micro3_p2e1 <- c42b_micro3_p2e1 %>% select(V4,V5,V6)
colnames(c42b_micro3_p2e1)[1] <- "V1"
colnames(c42b_micro3_p2e1)[2] <- "V2"
colnames(c42b_micro3_p2e1)[3] <- "V3"

c42b_micro_3_promoter2_with_ndr <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p2e1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter2_with_ndr <-rbind(c42b_micro_3_promoter2_with_ndr,some_count)
}

c42b_pn <- rbind(c42b_micro_3_promoter1_with_ndr, c42b_micro_3_promoter2_with_ndr)
c42b_pn <- unique(c42b_pn)

write.table(c42b_pn, "Promoter_in_loop_with_ndr.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")





c42b_micro3_p1e2 <- c42b_micro_3 %>% filter(Promoter_1 == 1 & H3K27me3_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_p2e1 <- c42b_micro_3 %>% filter(Promoter_2 == 1 & H3K27me3_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_promoter1_with_H3K27me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p1e2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter1_with_H3K27me3 <-rbind(c42b_micro_3_promoter1_with_H3K27me3,some_count)
}

c42b_micro3_p2e1 <- c42b_micro3_p2e1 %>% select(V4,V5,V6)
colnames(c42b_micro3_p2e1)[1] <- "V1"
colnames(c42b_micro3_p2e1)[2] <- "V2"
colnames(c42b_micro3_p2e1)[3] <- "V3"

c42b_micro_3_promoter2_with_H3K27me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p2e1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter2_with_H3K27me3 <-rbind(c42b_micro_3_promoter2_with_H3K27me3,some_count)
}

c42b_pe <- rbind(c42b_micro_3_promoter1_with_H3K27me3, c42b_micro_3_promoter2_with_H3K27me3)
c42b_pe <- unique(c42b_pe)

write.table(c42b_pe, "Promoter_in_loop_with_H3K27me3.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


c42b_micro3_p1h92 <- c42b_micro_3 %>% filter(Promoter_1 == 1 & H3K9me3_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_p2h91 <- c42b_micro_3 %>% filter(Promoter_2 == 1 & H3K9me3_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_promoter1_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p1h92 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter1_with_h3k9me3 <-rbind(c42b_micro_3_promoter1_with_h3k9me3,some_count)
}

c42b_micro3_p2h91 <- c42b_micro3_p2h91 %>% select(V4,V5,V6)
colnames(c42b_micro3_p2h91)[1] <- "V1"
colnames(c42b_micro3_p2h91)[2] <- "V2"
colnames(c42b_micro3_p2h91)[3] <- "V3"

c42b_micro_3_promoter2_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p2h91 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter2_with_h3k9me3 <-rbind(c42b_micro_3_promoter2_with_h3k9me3,some_count)
}

c42b_ph9 <- rbind(c42b_micro_3_promoter1_with_h3k9me3, c42b_micro_3_promoter2_with_h3k9me3)
c42b_ph9 <- unique(c42b_ph9)

write.table(c42b_ph9, "Promoter_in_loop_with_h3k9me3.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro3_p1x2 <- c42b_micro_3 %>% filter(Promoter_1 == 1 & None_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_p2x1 <- c42b_micro_3 %>% filter(Promoter_2 == 1 & None_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_promoter1_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p1x2 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter1_with_none <-rbind(c42b_micro_3_promoter1_with_none,some_count)
}

c42b_micro3_p2x1 <- c42b_micro3_p2x1 %>% select(V4,V5,V6)
colnames(c42b_micro3_p2x1)[1] <- "V1"
colnames(c42b_micro3_p2x1)[2] <- "V2"
colnames(c42b_micro3_p2x1)[3] <- "V3"

c42b_micro_3_promoter2_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_p2x1 %>% filter(V1 == i)
  chr_fithic <- h3k4me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_promoter2_with_none <-rbind(c42b_micro_3_promoter2_with_none,some_count)
}

c42b_px <- rbind(c42b_micro_3_promoter1_with_none, c42b_micro_3_promoter2_with_none)
c42b_px <- unique(c42b_px)

write.table(c42b_px, "Promoter_in_loop_with_none.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



c42b_micro3_e1p2 <- c42b_micro_3 %>% filter(Enhancer_1 == 1 & Promoter_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_e2p1 <- c42b_micro_3 %>% filter(Enhancer_2 == 1 & Promoter_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_enhancer1_with_promoter <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e1p2 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer1_with_promoter <-rbind(c42b_micro_3_enhancer1_with_promoter,some_count)
}

c42b_micro3_e2p1 <- c42b_micro3_e2p1 %>% select(V4,V5,V6)
colnames(c42b_micro3_e2p1)[1] <- "V1"
colnames(c42b_micro3_e2p1)[2] <- "V2"
colnames(c42b_micro3_e2p1)[3] <- "V3"

c42b_micro_3_enhancer2_with_promoter <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e2p1 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer2_with_promoter <-rbind(c42b_micro_3_enhancer2_with_promoter,some_count)
}

c42b_ep <- rbind(c42b_micro_3_enhancer1_with_promoter, c42b_micro_3_enhancer2_with_promoter)
c42b_ep <- unique(c42b_ep)

write.table(c42b_ep, "Enhancer_in_loop_with_promoter.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro3_e1i2 <- c42b_micro_3 %>% filter(Enhancer_1 == 1 & Insulator_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_e2i1 <- c42b_micro_3 %>% filter(Enhancer_2 == 1 & Insulator_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_enhancer1_with_insulator <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e1i2 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer1_with_insulator <-rbind(c42b_micro_3_enhancer1_with_insulator,some_count)
}

c42b_micro3_e2i1 <- c42b_micro3_e2i1 %>% select(V4,V5,V6)
colnames(c42b_micro3_e2i1)[1] <- "V1"
colnames(c42b_micro3_e2i1)[2] <- "V2"
colnames(c42b_micro3_e2i1)[3] <- "V3"

c42b_micro_3_enhancer2_with_insulator <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e2i1 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer2_with_insulator <-rbind(c42b_micro_3_enhancer2_with_insulator,some_count)
}

c42b_ei <- rbind(c42b_micro_3_enhancer1_with_insulator, c42b_micro_3_enhancer2_with_insulator)
c42b_ei <- unique(c42b_ei)

write.table(c42b_ei, "Enhancer_in_loop_with_insulator.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro3_e1n2 <- c42b_micro_3 %>% filter(Enhancer_1 == 1 & NDR_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_e2n1 <- c42b_micro_3 %>% filter(Enhancer_2 == 1 & NDR_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_enhancer1_with_ndr <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e1n2 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer1_with_ndr <-rbind(c42b_micro_3_enhancer1_with_ndr,some_count)
}

c42b_micro3_e2n1 <- c42b_micro3_e2n1 %>% select(V4,V5,V6)
colnames(c42b_micro3_e2n1)[1] <- "V1"
colnames(c42b_micro3_e2n1)[2] <- "V2"
colnames(c42b_micro3_e2n1)[3] <- "V3"

c42b_micro_3_enhancer2_with_ndr <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e2n1 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer2_with_ndr <-rbind(c42b_micro_3_enhancer2_with_ndr,some_count)
}

c42b_en <- rbind(c42b_micro_3_enhancer1_with_ndr, c42b_micro_3_enhancer2_with_ndr)
c42b_en <- unique(c42b_en)

write.table(c42b_en, "Enhancer_in_loop_with_ndr.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro3_e1h272 <- c42b_micro_3 %>% filter(Enhancer_1 == 1 & H3K27me3_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_e2h271 <- c42b_micro_3 %>% filter(Enhancer_2 == 1 & H3K27me3_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_enhancer1_with_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e1h272 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer1_with_h3k27me3 <-rbind(c42b_micro_3_enhancer1_with_h3k27me3,some_count)
}

c42b_micro3_e2h271 <- c42b_micro3_e2h271 %>% select(V4,V5,V6)
colnames(c42b_micro3_e2h271)[1] <- "V1"
colnames(c42b_micro3_e2h271)[2] <- "V2"
colnames(c42b_micro3_e2h271)[3] <- "V3"

c42b_micro_3_enhancer2_with_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e2h271 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer2_with_h3k27me3 <-rbind(c42b_micro_3_enhancer2_with_h3k27me3,some_count)
}

c42b_eh27 <- rbind(c42b_micro_3_enhancer1_with_h3k27me3, c42b_micro_3_enhancer2_with_h3k27me3)
c42b_eh27 <- unique(c42b_eh27)

write.table(c42b_eh27, "Enhancer_in_loop_with_h3k27me3.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro3_e1h92 <- c42b_micro_3 %>% filter(Enhancer_1 == 1 & H3K9me3_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_e2h91 <- c42b_micro_3 %>% filter(Enhancer_2 == 1 & H3K9me3_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_enhancer1_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e1h92 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer1_with_h3k9me3 <-rbind(c42b_micro_3_enhancer1_with_h3k9me3,some_count)
}

c42b_micro3_e2h91 <- c42b_micro3_e2h91 %>% select(V4,V5,V6)
colnames(c42b_micro3_e2h91)[1] <- "V1"
colnames(c42b_micro3_e2h91)[2] <- "V2"
colnames(c42b_micro3_e2h91)[3] <- "V3"

c42b_micro_3_enhancer2_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e2h91 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer2_with_h3k9me3 <-rbind(c42b_micro_3_enhancer2_with_h3k9me3,some_count)
}

c42b_eh9 <- rbind(c42b_micro_3_enhancer1_with_h3k9me3, c42b_micro_3_enhancer2_with_h3k9me3)
c42b_eh9 <- unique(c42b_eh9)

write.table(c42b_eh9, "Enhancer_in_loop_with_h3k9me3.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro3_e1x2 <- c42b_micro_3 %>% filter(Enhancer_1 == 1 & None_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro3_e2x1 <- c42b_micro_3 %>% filter(Enhancer_2 == 1 & None_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_enhancer1_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e1x2 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer1_with_none <-rbind(c42b_micro_3_enhancer1_with_none,some_count)
}

c42b_micro3_e2x1 <- c42b_micro3_e2x1 %>% select(V4,V5,V6)
colnames(c42b_micro3_e2x1)[1] <- "V1"
colnames(c42b_micro3_e2x1)[2] <- "V2"
colnames(c42b_micro3_e2x1)[3] <- "V3"

c42b_micro_3_enhancer2_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_e2x1 %>% filter(V1 == i)
  chr_fithic <- h3k27ac %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_enhancer2_with_none <-rbind(c42b_micro_3_enhancer2_with_none,some_count)
}

c42b_ex <- rbind(c42b_micro_3_enhancer1_with_none, c42b_micro_3_enhancer2_with_none)
c42b_ex <- unique(c42b_ex)

write.table(c42b_ex, "Enhancer_in_loop_with_none.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro3_i1p2 <- c42b_micro_3 %>% filter(Insulator_1 == 1 & Promoter_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_i2p1 <- c42b_micro_3 %>% filter(Insulator_2 == 1 & Promoter_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_insulator1_with_promoter <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_i1p2 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator1_with_promoter <-rbind(c42b_micro_3_insulator1_with_promoter,some_count)
}


c42b_micro_3_i2p1 <- c42b_micro_3_i2p1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_i2p1)[1] <- "V1"
colnames(c42b_micro_3_i2p1)[2] <- "V2"
colnames(c42b_micro_3_i2p1)[3] <- "V3"

c42b_micro_3_insulator2_with_promoter <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i2p1 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator2_with_promoter <-rbind(c42b_micro_3_insulator2_with_promoter,some_count)
}

c42b_ip <- rbind(c42b_micro_3_insulator1_with_promoter, c42b_micro_3_insulator2_with_promoter)
c42b_ip <- unique(c42b_ip)

write.table(c42b_ip, "Insulator_in_loop_with_promoter.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_i1e2 <- c42b_micro_3 %>% filter(Insulator_1 == 1 & Enhancer_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_i2e1 <- c42b_micro_3 %>% filter(Insulator_2 == 1 & Enhancer_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_insulator1_with_enhancer <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i1e2 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator1_with_enhancer <-rbind(c42b_micro_3_insulator1_with_enhancer,some_count)
}

c42b_micro_3_i2e1 <- c42b_micro_3_i2e1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_i2e1)[1] <- "V1"
colnames(c42b_micro_3_i2e1)[2] <- "V2"
colnames(c42b_micro_3_i2e1)[3] <- "V3"

c42b_micro_3_insulator2_with_enhancer <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i2e1 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator2_with_enhancer <-rbind(c42b_micro_3_insulator2_with_enhancer,some_count)
}

c42b_ie <- rbind(c42b_micro_3_insulator1_with_enhancer, c42b_micro_3_insulator2_with_enhancer)
c42b_ie <- unique(c42b_ie)

write.table(c42b_ie, "Insulator_in_loop_with_enhancer.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_i1n2 <- c42b_micro_3 %>% filter(Insulator_1 == 1 & NDR_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_i2n1 <- c42b_micro_3 %>% filter(Insulator_2 == 1 & NDR_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_insulator1_with_ndr <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i1n2 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator1_with_ndr <-rbind(c42b_micro_3_insulator1_with_ndr,some_count)
}

c42b_micro_3_i2n1 <- c42b_micro_3_i2n1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_i2n1)[1] <- "V1"
colnames(c42b_micro_3_i2n1)[2] <- "V2"
colnames(c42b_micro_3_i2n1)[3] <- "V3"

c42b_micro_3_insulator2_with_ndr <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i2n1 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator2_with_ndr <-rbind(c42b_micro_3_insulator2_with_ndr,some_count)
}

c42b_in <- rbind(c42b_micro_3_insulator1_with_ndr, c42b_micro_3_insulator2_with_ndr)
c42b_in <- unique(c42b_in)

write.table(c42b_in, "Insulator_in_loop_with_ndr.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_i1h272 <- c42b_micro_3 %>% filter(Insulator_1 == 1 & H3K27me3_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_i2h271 <- c42b_micro_3 %>% filter(Insulator_2 == 1 & H3K27me3_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_insulator1_with_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i1h272 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator1_with_h3k27me3 <-rbind(c42b_micro_3_insulator1_with_h3k27me3,some_count)
}

c42b_micro_3_i2h271 <- c42b_micro_3_i2h271 %>% select(V4,V5,V6)
colnames(c42b_micro_3_i2h271)[1] <- "V1"
colnames(c42b_micro_3_i2h271)[2] <- "V2"
colnames(c42b_micro_3_i2h271)[3] <- "V3"

c42b_micro_3_insulator2_with_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i2h271 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator2_with_h3k27me3 <-rbind(c42b_micro_3_insulator2_with_h3k27me3,some_count)
}

c42b_ih27 <- rbind(c42b_micro_3_insulator1_with_h3k27me3, c42b_micro_3_insulator2_with_h3k27me3)
c42b_ih27 <- unique(c42b_ih27)

write.table(c42b_ih27, "Insulator_in_loop_with_h3k27me3.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_i1h92 <- c42b_micro_3 %>% filter(Insulator_1 == 1 & H3K9me3_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_i2h91 <- c42b_micro_3 %>% filter(Insulator_2 == 1 & H3K9me3_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_insulator1_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i1h92 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator1_with_h3k9me3 <-rbind(c42b_micro_3_insulator1_with_h3k9me3,some_count)
}

c42b_micro_3_i2h91 <- c42b_micro_3_i2h91 %>% select(V4,V5,V6)
colnames(c42b_micro_3_i2h91)[1] <- "V1"
colnames(c42b_micro_3_i2h91)[2] <- "V2"
colnames(c42b_micro_3_i2h91)[3] <- "V3"

c42b_micro_3_insulator2_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i2h91 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator2_with_h3k9me3 <-rbind(c42b_micro_3_insulator2_with_h3k9me3,some_count)
}

c42b_ih9 <- rbind(c42b_micro_3_insulator1_with_h3k9me3, c42b_micro_3_insulator2_with_h3k9me3)
c42b_ih9 <- unique(c42b_ih9)

write.table(c42b_ih9, "Insulator_in_loop_with_h3k9me3.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_i1x2 <- c42b_micro_3 %>% filter(Insulator_1 == 1 & None_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_i2x1 <- c42b_micro_3 %>% filter(Insulator_2 == 1 & None_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_insulator1_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i1x2 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator1_with_none <-rbind(c42b_micro_3_insulator1_with_none,some_count)
}

c42b_micro_3_i2x1 <- c42b_micro_3_i2x1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_i2x1)[1] <- "V1"
colnames(c42b_micro_3_i2x1)[2] <- "V2"
colnames(c42b_micro_3_i2x1)[3] <- "V3"

c42b_micro_3_insulator2_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_i2x1 %>% filter(V1 == i)
  chr_fithic <- ctcf %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_insulator2_with_none <-rbind(c42b_micro_3_insulator2_with_none,some_count)
}

c42b_ix <- rbind(c42b_micro_3_insulator1_with_none, c42b_micro_3_insulator2_with_none)
c42b_ix <- unique(c42b_ix)

write.table(c42b_ix, "Insulator_in_loop_with_none.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro3_n1p2 <- c42b_micro_3 %>% filter(NDR_1 == 1 & Promoter_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2p1 <- c42b_micro_3 %>% filter(NDR_2 == 1 & Promoter_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_ndr1_with_promoter <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_n1p2 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr1_with_promoter <-rbind(c42b_micro_3_ndr1_with_promoter,some_count)
}

c42b_micro_3_n2p1 <- c42b_micro_3_n2p1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2p1)[1] <- "V1"
colnames(c42b_micro_3_n2p1)[2] <- "V2"
colnames(c42b_micro_3_n2p1)[3] <- "V3"

c42b_micro_3_ndr2_with_promoter <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2p1 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr2_with_promoter <-rbind(c42b_micro_3_ndr2_with_promoter,some_count)
}

c42b_np <- rbind(c42b_micro_3_ndr1_with_promoter, c42b_micro_3_ndr2_with_promoter)
c42b_np <- unique(c42b_np)

write.table(c42b_np, "NDR_in_loop_with_promoter.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1e2 <- c42b_micro_3 %>% filter(NDR_1 == 1 & Enhancer_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2e1 <- c42b_micro_3 %>% filter(NDR_2 == 1 & Enhancer_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_ndr1_with_enhancer <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1e2 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr1_with_enhancer <-rbind(c42b_micro_3_ndr1_with_enhancer,some_count)
}

c42b_micro_3_n2e1 <- c42b_micro_3_n2e1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2e1)[1] <- "V1"
colnames(c42b_micro_3_n2e1)[2] <- "V2"
colnames(c42b_micro_3_n2e1)[3] <- "V3"

c42b_micro_3_ndr2_with_enhancer <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2e1 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr2_with_enhancer <-rbind(c42b_micro_3_ndr2_with_enhancer,some_count)
}

c42b_ne <- rbind(c42b_micro_3_ndr1_with_enhancer, c42b_micro_3_ndr2_with_enhancer)
c42b_ne <- unique(c42b_ne)

write.table(c42b_ne, "NDR_in_loop_with_enhancer.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1i2 <- c42b_micro_3 %>% filter(NDR_1 == 1 & Insulator_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2i1 <- c42b_micro_3 %>% filter(NDR_2 == 1 & Insulator_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_ndr1_with_insulator <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1i2 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr1_with_insulator <-rbind(c42b_micro_3_ndr1_with_insulator,some_count)
}

c42b_micro_3_n2i1 <- c42b_micro_3_n2i1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2i1)[1] <- "V1"
colnames(c42b_micro_3_n2i1)[2] <- "V2"
colnames(c42b_micro_3_n2i1)[3] <- "V3"

c42b_micro_3_ndr2_with_insulator <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2i1 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr2_with_insulator <-rbind(c42b_micro_3_ndr2_with_insulator,some_count)
}

c42b_ni <- rbind(c42b_micro_3_ndr1_with_insulator, c42b_micro_3_ndr2_with_insulator)
c42b_ni <- unique(c42b_ni)

write.table(c42b_ni, "NDR_in_loop_with_insulator.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1h272 <- c42b_micro_3 %>% filter(NDR_1 == 1 & H3K27me3_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2h271 <- c42b_micro_3 %>% filter(NDR_2 == 1 & H3K27me3_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_ndr1_with_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1h272 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr1_with_h3k27me3 <-rbind(c42b_micro_3_ndr1_with_h3k27me3,some_count)
}

c42b_micro_3_n2h271 <- c42b_micro_3_n2h271 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2h271)[1] <- "V1"
colnames(c42b_micro_3_n2h271)[2] <- "V2"
colnames(c42b_micro_3_n2h271)[3] <- "V3"

c42b_micro_3_ndr2_with_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2h271 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr2_with_h3k27me3 <-rbind(c42b_micro_3_ndr2_with_h3k27me3,some_count)
}

c42b_nh27 <- rbind(c42b_micro_3_ndr1_with_h3k27me3, c42b_micro_3_ndr2_with_h3k27me3)
c42b_nh27 <- unique(c42b_nh27)

write.table(c42b_nh27, "NDR_in_loop_with_h3k27me3.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1h92 <- c42b_micro_3 %>% filter(NDR_1 == 1 & H3K9me3_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2h91 <- c42b_micro_3 %>% filter(NDR_2 == 1 & H3K9me3_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_ndr1_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1h92 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr1_with_h3k9me3 <-rbind(c42b_micro_3_ndr1_with_h3k9me3,some_count)
}

c42b_micro_3_n2h91 <- c42b_micro_3_n2h91 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2h91)[1] <- "V1"
colnames(c42b_micro_3_n2h91)[2] <- "V2"
colnames(c42b_micro_3_n2h91)[3] <- "V3"

c42b_micro_3_ndr2_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2h91 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr2_with_h3k9me3 <-rbind(c42b_micro_3_ndr2_with_h3k9me3,some_count)
}

c42b_nh9 <- rbind(c42b_micro_3_ndr1_with_h3k9me3, c42b_micro_3_ndr2_with_h3k9me3)
c42b_nh9 <- unique(c42b_nh9)

write.table(c42b_nh9, "NDR_in_loop_with_h3k9me3.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1x2 <- c42b_micro_3 %>% filter(NDR_1 == 1 & None_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2x1 <- c42b_micro_3 %>% filter(NDR_2 == 1 & None_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_ndr1_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1x2 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr1_with_none <-rbind(c42b_micro_3_ndr1_with_none,some_count)
}

c42b_micro_3_n2x1 <- c42b_micro_3_n2x1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2x1)[1] <- "V1"
colnames(c42b_micro_3_n2x1)[2] <- "V2"
colnames(c42b_micro_3_n2x1)[3] <- "V3"

c42b_micro_3_ndr2_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2x1 %>% filter(V1 == i)
  chr_fithic <- ndr %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_ndr2_with_none <-rbind(c42b_micro_3_ndr2_with_none,some_count)
}

c42b_nx <- rbind(c42b_micro_3_ndr1_with_none, c42b_micro_3_ndr2_with_none)
c42b_nx <- unique(c42b_nx)

write.table(c42b_nx, "NDR_in_loop_with_none.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro3_n1p2 <- c42b_micro_3 %>% filter(H3K27me3_1 == 1 & Promoter_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2p1 <- c42b_micro_3 %>% filter(H3K27me3_2 == 1 & Promoter_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k27me31_with_promoter <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_n1p2 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me31_with_promoter <-rbind(c42b_micro_3_h3k27me31_with_promoter,some_count)
}

c42b_micro_3_n2p1 <- c42b_micro_3_n2p1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2p1)[1] <- "V1"
colnames(c42b_micro_3_n2p1)[2] <- "V2"
colnames(c42b_micro_3_n2p1)[3] <- "V3"

c42b_micro_3_h3k27me32_with_promoter <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2p1 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me32_with_promoter <-rbind(c42b_micro_3_h3k27me32_with_promoter,some_count)
}

c42b_np <- rbind(c42b_micro_3_h3k27me31_with_promoter, c42b_micro_3_h3k27me32_with_promoter)
c42b_np <- unique(c42b_np)

write.table(c42b_np, "H3K27me3_in_loop_with_promoter.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1e2 <- c42b_micro_3 %>% filter(H3K27me3_1 == 1 & Enhancer_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2e1 <- c42b_micro_3 %>% filter(H3K27me3_2 == 1 & Enhancer_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k27me31_with_enhancer <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1e2 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me31_with_enhancer <-rbind(c42b_micro_3_h3k27me31_with_enhancer,some_count)
}

c42b_micro_3_n2e1 <- c42b_micro_3_n2e1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2e1)[1] <- "V1"
colnames(c42b_micro_3_n2e1)[2] <- "V2"
colnames(c42b_micro_3_n2e1)[3] <- "V3"

c42b_micro_3_h3k27me32_with_enhancer <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2e1 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me32_with_enhancer <-rbind(c42b_micro_3_h3k27me32_with_enhancer,some_count)
}

c42b_ne <- rbind(c42b_micro_3_h3k27me31_with_enhancer, c42b_micro_3_h3k27me32_with_enhancer)
c42b_ne <- unique(c42b_ne)

write.table(c42b_ne, "H3K27me3_in_loop_with_enhancer.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1i2 <- c42b_micro_3 %>% filter(H3K27me3_1 == 1 & Insulator_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2i1 <- c42b_micro_3 %>% filter(H3K27me3_2 == 1 & Insulator_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k27me31_with_insulator <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1i2 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me31_with_insulator <-rbind(c42b_micro_3_h3k27me31_with_insulator,some_count)
}

c42b_micro_3_n2i1 <- c42b_micro_3_n2i1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2i1)[1] <- "V1"
colnames(c42b_micro_3_n2i1)[2] <- "V2"
colnames(c42b_micro_3_n2i1)[3] <- "V3"

c42b_micro_3_h3k27me32_with_insulator <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2i1 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me32_with_insulator <-rbind(c42b_micro_3_h3k27me32_with_insulator,some_count)
}

c42b_ni <- rbind(c42b_micro_3_h3k27me31_with_insulator, c42b_micro_3_h3k27me32_with_insulator)
c42b_ni <- unique(c42b_ni)

write.table(c42b_ni, "H3K27me3_in_loop_with_insulator.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1h272 <- c42b_micro_3 %>% filter(H3K27me3_1 == 1 & NDR_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2h271 <- c42b_micro_3 %>% filter(H3K27me3_2 == 1 & NDR_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k27me31_with_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1h272 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me31_with_h3k27me3 <-rbind(c42b_micro_3_h3k27me31_with_h3k27me3,some_count)
}

c42b_micro_3_n2h271 <- c42b_micro_3_n2h271 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2h271)[1] <- "V1"
colnames(c42b_micro_3_n2h271)[2] <- "V2"
colnames(c42b_micro_3_n2h271)[3] <- "V3"

c42b_micro_3_h3k27me32_with_h3k27me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2h271 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me32_with_h3k27me3 <-rbind(c42b_micro_3_h3k27me32_with_h3k27me3,some_count)
}

c42b_nh27 <- rbind(c42b_micro_3_h3k27me31_with_h3k27me3, c42b_micro_3_h3k27me32_with_h3k27me3)
c42b_nh27 <- unique(c42b_nh27)

write.table(c42b_nh27, "H3K27me3_in_loop_with_NDR.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1h92 <- c42b_micro_3 %>% filter(H3K27me3_1 == 1 & H3K9me3_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2h91 <- c42b_micro_3 %>% filter(H3K27me3_2 == 1 & H3K9me3_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k27me31_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1h92 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me31_with_h3k9me3 <-rbind(c42b_micro_3_h3k27me31_with_h3k9me3,some_count)
}

c42b_micro_3_n2h91 <- c42b_micro_3_n2h91 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2h91)[1] <- "V1"
colnames(c42b_micro_3_n2h91)[2] <- "V2"
colnames(c42b_micro_3_n2h91)[3] <- "V3"

c42b_micro_3_h3k27me32_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2h91 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me32_with_h3k9me3 <-rbind(c42b_micro_3_h3k27me32_with_h3k9me3,some_count)
}

c42b_nh9 <- rbind(c42b_micro_3_h3k27me31_with_h3k9me3, c42b_micro_3_h3k27me32_with_h3k9me3)
c42b_nh9 <- unique(c42b_nh9)

write.table(c42b_nh9, "H3K27me3_in_loop_with_h3k9me3.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1x2 <- c42b_micro_3 %>% filter(H3K27me3_1 == 1 & None_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2x1 <- c42b_micro_3 %>% filter(H3K27me3_2 == 1 & None_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k27me31_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1x2 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me31_with_none <-rbind(c42b_micro_3_h3k27me31_with_none,some_count)
}

c42b_micro_3_n2x1 <- c42b_micro_3_n2x1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2x1)[1] <- "V1"
colnames(c42b_micro_3_n2x1)[2] <- "V2"
colnames(c42b_micro_3_n2x1)[3] <- "V3"

c42b_micro_3_h3k27me32_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2x1 %>% filter(V1 == i)
  chr_fithic <- h3k27me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k27me32_with_none <-rbind(c42b_micro_3_h3k27me32_with_none,some_count)
}

c42b_nx <- rbind(c42b_micro_3_h3k27me31_with_none, c42b_micro_3_h3k27me32_with_none)
c42b_nx <- unique(c42b_nx)

write.table(c42b_nx, "H3K27me3_in_loop_with_none.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro3_n1p2 <- c42b_micro_3 %>% filter(H3K9me3_1 == 1 & Promoter_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2p1 <- c42b_micro_3 %>% filter(H3K9me3_2 == 1 & Promoter_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k9me31_with_promoter <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro3_n1p2 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me31_with_promoter <-rbind(c42b_micro_3_h3k9me31_with_promoter,some_count)
}

c42b_micro_3_n2p1 <- c42b_micro_3_n2p1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2p1)[1] <- "V1"
colnames(c42b_micro_3_n2p1)[2] <- "V2"
colnames(c42b_micro_3_n2p1)[3] <- "V3"

c42b_micro_3_h3k9me32_with_promoter <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2p1 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me32_with_promoter <-rbind(c42b_micro_3_h3k9me32_with_promoter,some_count)
}

c42b_np <- rbind(c42b_micro_3_h3k9me31_with_promoter, c42b_micro_3_h3k9me32_with_promoter)
c42b_np <- unique(c42b_np)

write.table(c42b_np, "H3K9me3_in_loop_with_promoter.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1e2 <- c42b_micro_3 %>% filter(H3K9me3_1 == 1 & Enhancer_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2e1 <- c42b_micro_3 %>% filter(H3K9me3_2 == 1 & Enhancer_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k9me31_with_enhancer <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1e2 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me31_with_enhancer <-rbind(c42b_micro_3_h3k9me31_with_enhancer,some_count)
}

c42b_micro_3_n2e1 <- c42b_micro_3_n2e1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2e1)[1] <- "V1"
colnames(c42b_micro_3_n2e1)[2] <- "V2"
colnames(c42b_micro_3_n2e1)[3] <- "V3"

c42b_micro_3_h3k9me32_with_enhancer <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2e1 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me32_with_enhancer <-rbind(c42b_micro_3_h3k9me32_with_enhancer,some_count)
}

c42b_ne <- rbind(c42b_micro_3_h3k9me31_with_enhancer, c42b_micro_3_h3k9me32_with_enhancer)
c42b_ne <- unique(c42b_ne)

write.table(c42b_ne, "H3K9me3_in_loop_with_enhancer.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1i2 <- c42b_micro_3 %>% filter(H3K9me3_1 == 1 & Insulator_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2i1 <- c42b_micro_3 %>% filter(H3K9me3_2 == 1 & Insulator_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k9me31_with_insulator <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1i2 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me31_with_insulator <-rbind(c42b_micro_3_h3k9me31_with_insulator,some_count)
}

c42b_micro_3_n2i1 <- c42b_micro_3_n2i1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2i1)[1] <- "V1"
colnames(c42b_micro_3_n2i1)[2] <- "V2"
colnames(c42b_micro_3_n2i1)[3] <- "V3"

c42b_micro_3_h3k9me32_with_insulator <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2i1 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me32_with_insulator <-rbind(c42b_micro_3_h3k9me32_with_insulator,some_count)
}

c42b_ni <- rbind(c42b_micro_3_h3k9me31_with_insulator, c42b_micro_3_h3k9me32_with_insulator)
c42b_ni <- unique(c42b_ni)

write.table(c42b_ni, "H3K9me3_in_loop_with_insulator.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1h272 <- c42b_micro_3 %>% filter(H3K9me3_1 == 1 & NDR_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2h271 <- c42b_micro_3 %>% filter(H3K9me3_2 == 1 & NDR_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k9me31_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1h272 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me31_with_h3k9me3 <-rbind(c42b_micro_3_h3k9me31_with_h3k9me3,some_count)
}

c42b_micro_3_n2h271 <- c42b_micro_3_n2h271 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2h271)[1] <- "V1"
colnames(c42b_micro_3_n2h271)[2] <- "V2"
colnames(c42b_micro_3_n2h271)[3] <- "V3"

c42b_micro_3_h3k9me32_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2h271 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me32_with_h3k9me3 <-rbind(c42b_micro_3_h3k9me32_with_h3k9me3,some_count)
}

c42b_nh27 <- rbind(c42b_micro_3_h3k9me31_with_h3k9me3, c42b_micro_3_h3k9me32_with_h3k9me3)
c42b_nh27 <- unique(c42b_nh27)

write.table(c42b_nh27, "H3K9me3_in_loop_with_NDR.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1h92 <- c42b_micro_3 %>% filter(H3K9me3_1 == 1 & H3K27me3_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2h91 <- c42b_micro_3 %>% filter(H3K9me3_2 == 1 & H3K27me3_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k9me31_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1h92 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me31_with_h3k9me3 <-rbind(c42b_micro_3_h3k9me31_with_h3k9me3,some_count)
}

c42b_micro_3_n2h91 <- c42b_micro_3_n2h91 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2h91)[1] <- "V1"
colnames(c42b_micro_3_n2h91)[2] <- "V2"
colnames(c42b_micro_3_n2h91)[3] <- "V3"

c42b_micro_3_h3k9me32_with_h3k9me3 <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2h91 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me32_with_h3k9me3 <-rbind(c42b_micro_3_h3k9me32_with_h3k9me3,some_count)
}

c42b_nh9 <- rbind(c42b_micro_3_h3k9me31_with_h3k9me3, c42b_micro_3_h3k9me32_with_h3k9me3)
c42b_nh9 <- unique(c42b_nh9)

write.table(c42b_nh9, "H3K9me3_in_loop_with_h3k27me3.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

c42b_micro_3_n1x2 <- c42b_micro_3 %>% filter(H3K9me3_1 == 1 & None_2 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)
c42b_micro_3_n2x1 <- c42b_micro_3 %>% filter(H3K9me3_2 == 1 & None_1 == 1) %>% select(V1,V2,V3,V4,V5,V6,V7,V8)

c42b_micro_3_h3k9me31_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n1x2 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me31_with_none <-rbind(c42b_micro_3_h3k9me31_with_none,some_count)
}

c42b_micro_3_n2x1 <- c42b_micro_3_n2x1 %>% select(V4,V5,V6)
colnames(c42b_micro_3_n2x1)[1] <- "V1"
colnames(c42b_micro_3_n2x1)[2] <- "V2"
colnames(c42b_micro_3_n2x1)[3] <- "V3"

c42b_micro_3_h3k9me32_with_none <- data.frame()
for (i in x){
  chr_mustache <- c42b_micro_3_n2x1 %>% filter(V1 == i)
  chr_fithic <- h3k9me3 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic,chr_mustache, max_dist= 10000, by = c("V2","V3"))
  c42b_micro_3_h3k9me32_with_none <-rbind(c42b_micro_3_h3k9me32_with_none,some_count)
}

c42b_nx <- rbind(c42b_micro_3_h3k9me31_with_none, c42b_micro_3_h3k9me32_with_none)
c42b_nx <- unique(c42b_nx)

write.table(c42b_nx, "H3K9me3_in_loop_with_none.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


