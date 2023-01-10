library(plyr)
library(dplyr)
library(fuzzyjoin)

#This script is similar to loop category script, but filter loops based on the loop information they have, and save them based on the loop categories so we can compare the loop information.

loop_reg_analysis <- function(
loop_data,
label,
h3k27ac_bed = '/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/C42B_H3K27ac_no_tss.bed',
h3k4me3_bed = '/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/C42B_TSS_FPKM_greater_than_0.5.tsv',
CTCF_bed = '/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/C42B_CTCF_no_tss_no_enhancer.bed',
NDR_bed = '/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/C42B_NDR_without_features.bed',
h3k27me3_bed = '/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/H3K27me3_not_in_enhancer_insulator_promoter_ndr_region.bed',
h3k9me3_bed = '/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/H3K9me3_not_in_enhancer_insulator_promoter_ndr_h3k27me3_region.bed',
output_path
){

loop_data$loopname <- paste(label,rownames(loop_data),sep = "_")

h3k27ac <- read.table(h3k27ac_bed)
h3k4me3 <- read.table(h3k4me3_bed)
ctcf <- read.table(CTCF_bed)
ndr <- read.table(NDR_bed)
h3k27me3 <- read.table(h3k27me3_bed)
h3k9me3 <- read.table(h3k9me3_bed)

loop_data$V1 <- gsub("chr","",loop_data$V1)
loop_data$V4 <- gsub("chr","",loop_data$V4)
ctcf$V1 <- gsub("chr","",ctcf$V1)
h3k27ac$V1 <- gsub("chr","",h3k27ac$V1)
h3k4me3$V1 <- gsub("chr","",h3k4me3$V1)
ndr$V1 <- gsub("chr","",ndr$V1)
h3k27me3$V1 <- gsub("chr","",h3k27me3$V1)
h3k9me3$V1 <- gsub("chr","",h3k9me3$V1)

loop_data <- loop_data %>% mutate(
                                        Enhancer_1 = 0, Enhancer_2 = 0,
                                        Promoter_1 = 0, Promoter_2 = 0,
                                        Insulator_1 = 0, Insulator_2 = 0,
                                        NDR_1 = 0, NDR_2 = 0,
                                        h3k27me3_1 = 0, h3k27me3_2 = 0,
                                        h3k9me3_1 = 0, h3k9me3_2 = 0,
                                        None_1 = 1, None_2 = 1)

x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")

anchor1 <- loop_data %>% select(V1,V2,V3,loopname)
anchor2 <- loop_data %>% select(V4,V5,V6,loopname)
colnames(anchor2)[1] <- "V1"
colnames(anchor2)[2] <- "V2"
colnames(anchor2)[3] <- "V3"

#For promoter
anchor1_h3k4me3 <- data.frame()
anchor2_h3k4me3 <- data.frame()
non_anchor1_h3k4me3 <- data.frame()
non_anchor2_h3k4me3 <- data.frame()
for (i in x){
        chr_fithic <- h3k4me3 %>% filter(V1 == i)
        #For anchor1
  chr_mustache <- anchor1 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  anchor1_h3k4me3 <-rbind(anchor1_h3k4me3,some_count)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_anchor1_h3k4me3 <- rbind(non_anchor1_h3k4me3,some_count)

  #For anchor2
  chr_mustache <- anchor2 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  anchor2_h3k4me3 <-rbind(anchor2_h3k4me3,some_count)
  some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
  non_anchor2_h3k4me3 <- rbind(non_anchor2_h3k4me3,some_count)
        }

#integrate loop_data regarding promoter section
loop_data$Promoter_1[loop_data$loopname  %in% anchor1_h3k4me3$loopname] <- 1
loop_data$Promoter_2[loop_data$loopname  %in% anchor2_h3k4me3$loopname] <- 1
loop_data$None_1[loop_data$loopname  %in% anchor1_h3k4me3$loopname] <- 0
loop_data$None_2[loop_data$loopname  %in% anchor2_h3k4me3$loopname] <- 0


#For enhancer
anchor1_h3k27ac <- data.frame()
anchor2_h3k27ac <- data.frame()
non_anchor1_h3k27ac <- data.frame()
non_anchor2_h3k27ac <- data.frame()

for (i in x){
        chr_fithic <- h3k27ac %>% filter(V1 == i)
        #For anchor1
        chr_mustache <- non_anchor1_h3k4me3 %>% filter(V1 == i)
        some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        anchor1_h3k27ac <-rbind(anchor1_h3k27ac,some_count)
        some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        non_anchor1_h3k27ac <- rbind(non_anchor1_h3k27ac,some_count)

        #For anchor2
        chr_mustache <- non_anchor2_h3k4me3 %>% filter(V1 == i)
        some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        anchor2_h3k27ac <-rbind(anchor2_h3k27ac,some_count)
        some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        non_anchor2_h3k27ac <- rbind(non_anchor2_h3k27ac,some_count)
        }

#integrate loop_data regarding enhancer section
loop_data$Enhancer_1[loop_data$loopname  %in% anchor1_h3k27ac$loopname] <- 1
loop_data$Enhancer_2[loop_data$loopname  %in% anchor2_h3k27ac$loopname] <- 1
loop_data$None_1[loop_data$loopname  %in% anchor1_h3k27ac$loopname] <- 0
loop_data$None_2[loop_data$loopname  %in% anchor2_h3k27ac$loopname] <- 0

#For Insulator
anchor1_ctcf <- data.frame()
anchor2_ctcf <- data.frame()
non_anchor1_ctcf <- data.frame()
non_anchor2_ctcf <- data.frame()

for (i in x){
        chr_fithic <- ctcf %>% filter(V1 == i)
        #For anchor1
        chr_mustache <- non_anchor1_h3k27ac %>% filter(V1 == i)
        some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        anchor1_ctcf <-rbind(anchor1_ctcf,some_count)
        some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        non_anchor1_ctcf <- rbind(non_anchor1_ctcf,some_count)

        #For anchor2
        chr_mustache <- non_anchor2_h3k27ac %>% filter(V1 == i)
        some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        anchor2_ctcf <-rbind(anchor2_ctcf,some_count)
        some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        non_anchor2_ctcf <- rbind(non_anchor2_ctcf,some_count)
        }

#integrate loop_data regarding enhancer section
loop_data$Insulator_1[loop_data$loopname  %in% anchor1_ctcf$loopname] <- 1
loop_data$Insulator_2[loop_data$loopname  %in% anchor2_ctcf$loopname] <- 1
loop_data$None_1[loop_data$loopname  %in% anchor1_ctcf$loopname] <- 0
loop_data$None_2[loop_data$loopname  %in% anchor2_ctcf$loopname] <- 0

#For NDR
anchor1_ndr <- data.frame()
anchor2_ndr <- data.frame()
non_anchor1_ndr <- data.frame()
non_anchor2_ndr <- data.frame()

for (i in x){
        chr_fithic <- ndr %>% filter(V1 == i)
        #For anchor1
        chr_mustache <- non_anchor1_ctcf %>% filter(V1 == i)
        some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        anchor1_ndr <-rbind(anchor1_ndr,some_count)
        some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        non_anchor1_ndr <- rbind(non_anchor1_ndr,some_count)

        #For anchor2
        chr_mustache <- non_anchor2_ctcf %>% filter(V1 == i)
        some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        anchor2_ndr <-rbind(anchor2_ndr,some_count)
        some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        non_anchor2_ndr <- rbind(non_anchor2_ndr,some_count)
        }

#integrate loop_data regarding NDR section
loop_data$NDR_1[loop_data$loopname  %in% anchor1_ndr$loopname] <- 1
loop_data$NDR_2[loop_data$loopname  %in% anchor2_ndr$loopname] <- 1
loop_data$None_1[loop_data$loopname  %in% anchor1_ctcf$loopname] <- 0
loop_data$None_2[loop_data$loopname  %in% anchor2_ctcf$loopname] <- 0

#For repressed region - H3K27me3
anchor1_h3k27me3 <- data.frame()
anchor2_h3k27me3 <- data.frame()
non_anchor1_h3k27me3 <- data.frame()
non_anchor2_h3k27me3 <- data.frame()

for (i in x){
        chr_fithic <- h3k27me3 %>% filter(V1 == i)
        #For anchor1
        chr_mustache <- non_anchor1_ndr %>% filter(V1 == i)
        some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        anchor1_h3k27me3 <-rbind(anchor1_h3k27me3,some_count)
        some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        non_anchor1_h3k27me3 <- rbind(non_anchor1_h3k27me3,some_count)

        #For anchor2
        chr_mustache <- non_anchor2_ndr %>% filter(V1 == i)
        some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        anchor2_h3k27me3 <-rbind(anchor2_h3k27me3,some_count)
        some_count <- difference_anti_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        non_anchor2_h3k27me3 <- rbind(non_anchor2_h3k27me3,some_count)
        }

#integrate loop_data regarding repressed section
loop_data$H3K27me3_1[loop_data$loopname  %in% anchor1_h3k27me3$loopname] <- 1
loop_data$H3K27me3_2[loop_data$loopname  %in% anchor2_h3k27me3$loopname] <- 1
loop_data$None_1[loop_data$loopname  %in% anchor1_h3k27me3$loopname] <- 0
loop_data$None_2[loop_data$loopname  %in% anchor2_h3k27me3$loopname] <- 0

#For heterochromatin region - H3K9me3
anchor1_h3k9me3 <- data.frame()
anchor2_h3k9me3 <- data.frame()

for (i in x){
        chr_fithic <- h3k9me3 %>% filter(V1 == i)
        #For anchor1
        chr_mustache <- non_anchor1_h3k27me3 %>% filter(V1 == i)
        some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        anchor1_h3k9me3 <-rbind(anchor1_h3k9me3,some_count)

        #For anchor2
        chr_mustache <- non_anchor2_h3k27me3 %>% filter(V1 == i)
        some_count <- difference_semi_join(chr_mustache,chr_fithic, max_dist= 10000, by = c("V2","V3"))
        anchor2_h3k9me3 <-rbind(anchor2_h3k9me3,some_count)
        }

#integrate loop_data regarding heterochromatin section
loop_data$H3K9me3_1[loop_data$loopname  %in% anchor1_h3k9me3$loopname] <- 1
loop_data$H3K9me3_2[loop_data$loopname  %in% anchor2_h3k9me3$loopname] <- 1
loop_data$None_1[loop_data$loopname  %in% anchor1_h3k27me3$loopname] <- 0
loop_data$None_2[loop_data$loopname  %in% anchor2_h3k27me3$loopname] <- 0

#pick up the regulatory element interaction categories
loop_data_pp <- loop_data %>% filter(Promoter_1 == 1, Promoter_2 ==1)
loop_data_pe1 <- loop_data %>% filter(Promoter_1 == 1, Enhancer_2 ==1)
loop_data_pe2 <- loop_data %>% filter(Promoter_2 == 1, Enhancer_1 ==1)
loop_data_pe <- rbind(loop_data_pe1,loop_data_pe2)

loop_data_pi1 <- loop_data %>% filter(Promoter_1 == 1, Insulator_2 ==1)
loop_data_pi2 <- loop_data %>% filter(Promoter_2 == 1, Insulator_1 ==1)
loop_data_pi <- rbind(loop_data_pi1,loop_data_pi2)

loop_data_pn1 <- loop_data %>% filter(Promoter_1 == 1, NDR_2 ==1)
loop_data_pn2 <- loop_data %>% filter(Promoter_2 == 1, NDR_1 ==1)
loop_data_pn <- rbind(loop_data_pn1,loop_data_pn2)

loop_data_ph271 <- loop_data %>% filter(Promoter_1 == 1, H3K27me3_2 ==1)
loop_data_ph272 <- loop_data %>% filter(Promoter_2 == 1, H3K27me3_1 ==1)
loop_data_ph27 <- rbind(loop_data_ph271,loop_data_ph272)

loop_data_ph91 <- loop_data %>% filter(Promoter_1 == 1, H3K9me3_2 ==1)
loop_data_ph92 <- loop_data %>% filter(Promoter_2 == 1, H3K9me3_1 ==1)
loop_data_ph9 <- rbind(loop_data_ph91,loop_data_ph92)

loop_data_px1 <- loop_data %>% filter(Promoter_1 == 1, None_2 ==1)
loop_data_px2 <- loop_data %>% filter(Promoter_2 == 1, None_1 ==1)
loop_data_px <- rbind(loop_data_px1,loop_data_px2)

loop_data_ee <- loop_data %>% filter(Enhancer_1 == 1, Enhancer_2 ==1)

loop_data_ei1 <- loop_data %>% filter(Enhancer_1 == 1, Insulator_2 ==1)
loop_data_ei2 <- loop_data %>% filter(Enhancer_2 == 1, Insulator_1 ==1)
loop_data_ei <- rbind(loop_data_ei1,loop_data_ei2)

loop_data_en1 <- loop_data %>% filter(Enhancer_1 == 1, NDR_2 ==1)
loop_data_en2 <- loop_data %>% filter(Enhancer_2 == 1, NDR_1 ==1)
loop_data_en <- rbind(loop_data_en1,loop_data_en2)

loop_data_eh271 <- loop_data %>% filter(Enhancer_1 == 1, H3K27me3_2 ==1)
loop_data_eh272 <- loop_data %>% filter(Enhancer_2 == 1, H3K27me3_1 ==1)
loop_data_eh27 <- rbind(loop_data_eh271,loop_data_eh272)

loop_data_eh91 <- loop_data %>% filter(Enhancer_1 == 1, H3K9me3_2 ==1)
loop_data_eh92 <- loop_data %>% filter(Enhancer_2 == 1, H3K9me3_1 ==1)
loop_data_eh9 <- rbind(loop_data_eh91,loop_data_eh92)

loop_data_ex1 <- loop_data %>% filter(Enhancer_1 == 1, None_2 ==1)
loop_data_ex2 <- loop_data %>% filter(Enhancer_2 == 1, None_1 ==1)
loop_data_ex <- rbind(loop_data_ex1,loop_data_ex2)

loop_data_ii <- loop_data %>% filter(Insulator_1 == 1, Insulator_2 ==1)

loop_data_in1 <- loop_data %>% filter(Insulator_1 == 1, NDR_2 ==1)
loop_data_in2 <- loop_data %>% filter(Insulator_2 == 1, NDR_1 ==1)
loop_data_in <- rbind(loop_data_in1,loop_data_in2)

loop_data_ih271 <- loop_data %>% filter(Insulator_1 == 1, H3K27me3_2 ==1)
loop_data_ih272 <- loop_data %>% filter(Insulator_2 == 1, H3K27me3_1 ==1)
loop_data_ih27 <- rbind(loop_data_ih271,loop_data_ih272)

loop_data_ih91 <- loop_data %>% filter(Insulator_1 == 1, H3K9me3_2 ==1)
loop_data_ih92 <- loop_data %>% filter(Insulator_2 == 1, H3K9me3_1 ==1)
loop_data_ih9 <- rbind(loop_data_ih91,loop_data_ih92)

loop_data_ix1 <- loop_data %>% filter(Insulator_1 == 1, None_2 ==1)
loop_data_ix2 <- loop_data %>% filter(Insulator_2 == 1, None_1 ==1)
loop_data_ix <- rbind(loop_data_ix1,loop_data_ix2)

loop_data_nn <- loop_data %>% filter(NDR_1 == 1, NDR_2 ==1)


loop_data_nh271 <- loop_data %>% filter(NDR_1 == 1, H3K27me3_2 ==1)
loop_data_nh272 <- loop_data %>% filter(NDR_2 == 1, H3K27me3_1 ==1)
loop_data_nh27 <- rbind(loop_data_nh271,loop_data_nh272)

loop_data_nh91 <- loop_data %>% filter(NDR_1 == 1, H3K9me3_2 ==1)
loop_data_nh92 <- loop_data %>% filter(NDR_2 == 1, H3K9me3_1 ==1)
loop_data_nh9 <- rbind(loop_data_nh91,loop_data_nh92)

loop_data_nx1 <- loop_data %>% filter(NDR_1 == 1, None_2 ==1)
loop_data_nx2 <- loop_data %>% filter(NDR_2 == 1, None_1 ==1)
loop_data_nx <- rbind(loop_data_nx1,loop_data_nx2)

loop_data_h27h27 <- loop_data %>% filter(H3K27me3_1 == 1, H3K27me3_2 ==1)

loop_data_h27h91 <- loop_data %>% filter(H3K27me3_1 == 1, H3K9me3_2 ==1)
loop_data_h27h92 <- loop_data %>% filter(H3K27me3_2 == 1, H3K9me3_1 ==1)
loop_data_h27h9 <- rbind(loop_data_h27h91,loop_data_h27h92)

loop_data_h27x1 <- loop_data %>% filter(H3K27me3_1 == 1, None_2 ==1)
loop_data_h27x2 <- loop_data %>% filter(H3K27me3_2 == 1, None_1 ==1)
loop_data_h27x <- rbind(loop_data_h27x1,loop_data_h27x2)

loop_data_h9h9 <- loop_data %>% filter(H3K9me3_1 == 1, H3K9me3_2 ==1)

loop_data_h9x1 <- loop_data %>% filter(H3K9me3_1 == 1, None_2 ==1)
loop_data_h9x2 <- loop_data %>% filter(H3K9me3_2 == 1, None_1 ==1)
loop_data_h9x <- rbind(loop_data_h9x1,loop_data_h9x2)

loop_data_xx <- loop_data %>% filter(None_1 == 1, None_2 ==1)

write.table(loop_data_pp, file = paste0(output_path,'/',label,'_Promoter_Promoter_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_pe, file = paste0(output_path,'/',label,'_Promoter_Enhancer_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_pi, file = paste0(output_path,'/',label,'_Promoter_Insulator_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_pn, file = paste0(output_path,'/',label,'_Promoter_NDR_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_ph27, file = paste0(output_path,'/',label,'_Promoter_H3K27me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_ph9, file = paste0(output_path,'/',label,'_Promoter_H3K9me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_px, file = paste0(output_path,'/',label,'_Promoter_None_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(loop_data_ee, file = paste0(output_path,'/',label,'_Enhancer_Enhancer_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_ei, file = paste0(output_path,'/',label,'_Enhancer_Insulator_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_en, file = paste0(output_path,'/',label,'_Enhancer_NDR_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_eh27, file = paste0(output_path,'/',label,'_Enhancer_H3K27me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_eh9, file = paste0(output_path,'/',label,'_Enhancer_H3K9me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_ex, file = paste0(output_path,'/',label,'_Enhancer_None_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(loop_data_ii, file = paste0(output_path,'/',label,'_Insulator_Insulator_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_in, file = paste0(output_path,'/',label,'_Insulator_NDR_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_ih27, file = paste0(output_path,'/',label,'_Insulator_H3K27me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_ih9, file = paste0(output_path,'/',label,'_Insulator_H3K9me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_ix, file = paste0(output_path,'/',label,'_Insulator_None_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(loop_data_nn, file = paste0(output_path,'/',label,'_NDR_NDR_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_nh27, file = paste0(output_path,'/',label,'_NDR_H3K27me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_nh9, file = paste0(output_path,'/',label,'_NDR_H3K9me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_nx, file = paste0(output_path,'/',label,'_NDR_None_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(loop_data_h27h27, file = paste0(output_path,'/',label,'_H3K27me3_H3K27me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_h27h9, file = paste0(output_path,'/',label,'_H3K27me3_H3K9me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_h27x, file = paste0(output_path,'/',label,'_H3K27me3_None_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(loop_data_h9h9, file = paste0(output_path,'/',label,'_H3K9me3_H3K9me3_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(loop_data_h9x, file = paste0(output_path,'/',label,'_H3K9me3_None_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

write.table(loop_data_xx, file = paste0(output_path,'/',label,'_None_None_Loop.tsv'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

summarytable <- data.frame(
                        Reg_Categories = c('Total',
                        'Promoter_Promoter', 'Promoter_Enhancer', 'Promoter_Insulator', 'Promoter_NDR', 'Promoter_H3K27me3', 'Promoter_H3K9me3', 'Promoter_None',
                        'Enhancer_Enhancer', 'Enhancer_Insulator', 'Enhancer_NDR', 'Enhancer_H3K27me3', 'Enhancer_H3K9me3' , 'Enhancer_None',
                        'Insulator_Insulator', 'Insulator_NDR', 'Insulator_H3K27me3', 'Insulator_H3K9me3', 'Insulator_None',
                        'NDR_NDR', 'NDR_H3K27me3', 'NDR_H3K9me3', 'NDR_None',
                        'H3K27me3_H3K27me3', 'H3K27me3_H3K9me3', 'H3K27me3_None',
                        'H3K9me3_H3K9me3', 'H3K9me3_None',
                        'None_None'),
                        Number = c(nrow(loop_data),
                        nrow(loop_data_pp), nrow(loop_data_pe), nrow(loop_data_pi), nrow(loop_data_pn), nrow(loop_data_ph27), nrow(loop_data_ph9), nrow(loop_data_px),
                        nrow(loop_data_ee), nrow(loop_data_ei), nrow(loop_data_en), nrow(loop_data_eh27), nrow(loop_data_eh9), nrow(loop_data_ex),
                        nrow(loop_data_ii), nrow(loop_data_in), nrow(loop_data_ih27), nrow(loop_data_ih9), nrow(loop_data_ix),
                        nrow(loop_data_nn), nrow(loop_data_nh27), nrow(loop_data_nh9), nrow(loop_data_nx),
                        nrow(loop_data_h27h27), nrow(loop_data_h27h9), nrow(loop_data_h27x),
                        nrow(loop_data_h9h9), nrow(loop_data_h9x),
                        nrow(loop_data_xx)
                          )
                        ) #end data.frame
write.table(summarytable, file = paste0(output_path,'/',label,'_RegSummary.txt'), row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
}
