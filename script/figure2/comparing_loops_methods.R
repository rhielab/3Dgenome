library(dplyr)
library(plyr)
library(fuzzyjoin)
library(ggplot2)

#This script compare and identify the shared loops,and store that information as the matrices represented as 0 or 1 to show whether it shares any loop with the another methods.


#Load the loops
c42b_hic_5kb <-  read.table("HiC_Mustache-5kb-all-loop.tsv")
c42b_micro_1_5kb <- read.table("1_Billion_Mustache-5kb-all-loop.tsv")
c42b_micro_3_5kb <- read.table("3_billion_low_Mustache-5kb-all-loop.tsv")
c42b_micro_2_5kb <- read.table("2_Billion_Mustache-5kb-all-loop.tsv")

#Add the loop names with numbers; this is used as an identifier later on
c42b_micro_2_5kb$loopname <- paste("c42b_2_billion",rownames(c42b_micro_2_5kb),sep = "_")
c42b_hic_5kb$loopname <- paste("c42b_2_billion",rownames(c42b_hic_5kb),sep = "_")
c42b_micro_1_5kb$loopname <- paste("c42b_1_billion",rownames(c42b_micro_1_5kb),sep = "_")
c42b_micro_3_5kb$loopname <- paste("c42b_3_billion",rownames(c42b_micro_3_5kb),sep = "_")

#Add the score of 1 for its own methods, and 0 since it has yet to be identified as shared.
c42b_hic_5kb <- c42b_hic_5kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 1, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0)
c42b_micro_1_5kb <- c42b_micro_1_5kb %>% mutate(Micro_C_1_Billion = 1, Hi_C = 0, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0, Chicago = 0)
c42b_micro_3_5kb <- c42b_micro_3_5kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 0, Micro_C_2_Billion = 0,Micro_C_3_Billion = 1, Chicago = 0)
c42b_micro_2_5kb <- c42b_micro_2_5kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 0, Micro_C_2_Billion = 1,Micro_C_3_Billion = 0, Chicago = 0)

#This step is required to set all the methods equal to each other. Some have chr in their name, while some doesn't; by doing this, it makes so all of the methods have equivalent name to be compared later on.
c42b_micro_2_5kb$V1 <- gsub("chr","",c42b_micro_2_5kb$V1)
c42b_micro_2_5kb$V4 <- gsub("chr","",c42b_micro_2_5kb$V4)

c42b_chicago_5kb <-  read.table("PC-UNI1945.5kb.bedpe")
c42b_chicago_1kb <-  read.table("PC-UNI1945.1kb.bedpe")
c42b_chicago_2kb <-  read.table("PC-UNI1945.2kb.bedpe")
c42b_chicago_10kb <-  read.table("PC-UNI1945.10kb.bedpe")

c42b_chicago_5kb$V1 <- gsub("chr","",c42b_chicago_5kb$V1)
c42b_chicago_5kb$V4 <- gsub("chr","",c42b_chicago_5kb$V4)


c42b_chicago_1kb$V1 <- gsub("chr","",c42b_chicago_1kb$V1)
c42b_chicago_1kb$V4 <- gsub("chr","",c42b_chicago_1kb$V4)


c42b_chicago_2kb$V1 <- gsub("chr","",c42b_chicago_2kb$V1)
c42b_chicago_2kb$V4 <- gsub("chr","",c42b_chicago_2kb$V4)


c42b_chicago_10kb$V1 <- gsub("chr","",c42b_chicago_10kb$V1)
c42b_chicago_10kb$V4 <- gsub("chr","",c42b_chicago_10kb$V4)

c42b_chicago_5kb$loopname <- paste("c42b_chicago",rownames(c42b_chicago_5kb),sep = "_")
c42b_chicago_5kb <- c42b_chicago_5kb %>% mutate(Micro_C_1_Billion = 0, Chicago = 1, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0)
c42b_chicago_1kb$loopname <- paste("c42b_chicago",rownames(c42b_chicago_1kb),sep = "_")
c42b_chicago_2kb$loopname <- paste("c42b_chicago",rownames(c42b_chicago_2kb),sep = "_")
c42b_chicago_10kb$loopname <- paste("c42b_chicago",rownames(c42b_chicago_10kb),sep = "_")
c42b_chicago_1kb <- c42b_chicago_1kb %>% mutate(Micro_C_1_Billion = 0, Chicago = 1, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0)
c42b_chicago_2kb <- c42b_chicago_2kb %>% mutate(Micro_C_1_Billion = 0, Chicago = 1, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0)
c42b_chicago_10kb <- c42b_chicago_10kb %>% mutate(Micro_C_1_Billion = 0, Chicago = 1, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0)

#Chr to loop through
x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")

#Make empty data frame initially
micro_1_hic_5kb <- data.frame()
for (i in x){
#Filter based on chromosomes first, since difference semi join cannot interpret string(X/Y) and also will semi join chr1 and chr10
  chr_template <- c42b_micro_1_5kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_5kb %>% filter(V1 == i)
#Difference semi join extends to plus minus 10000 in the coordinate, so for example, if your coordinate is at 20000 25000 50000 55000 and there are three loops at 10000 15000 40000 45000, 30000 35000 60000 65000, and  50000 55000 20000 25000, first two will be recognized as similar loop, while the last one won't be, since it is more than 10000 away from the coordinates.
#I used 10000 for value here, to extend to two bins max. For others, you can also see this is the case, as I used 2000 for 1kb, 4000 for 2kb, etc.
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
#Add the shared loops, and bind them into the micro_1_hic_5kb
  micro_1_hic_5kb <-rbind(micro_1_hic_5kb,some_count)
}
#Check whether following loop name is in micro_1_hic_5kb, which is shared loop of Micro-C 1 Billion in Hi-C. If yes, change the number to 1, which means "loop is shared"
c42b_micro_1_5kb$Hi_C[c42b_micro_1_5kb$loopname %in% micro_1_hic_5kb$loopname] <- 1

micro_1_micro_2_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_1_micro_2_5kb <-rbind(micro_1_micro_2_5kb,some_count)
}
c42b_micro_1_5kb$Micro_C_2_Billion[c42b_micro_1_5kb$loopname %in% micro_1_micro_2_5kb$loopname] <- 1

micro_1_micro_3_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_1_micro_3_5kb <-rbind(micro_1_micro_3_5kb,some_count)
}
c42b_micro_1_5kb$Micro_C_3_Billion[c42b_micro_1_5kb$loopname %in% micro_1_micro_3_5kb$loopname] <- 1



micro_2_hic_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_5kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_2_hic_5kb <-rbind(micro_2_hic_5kb,some_count)
}
c42b_micro_2_5kb$Hi_C[c42b_micro_2_5kb$loopname %in% micro_2_hic_5kb$loopname] <- 1

micro_2_micro_1_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_2_micro_1_5kb <-rbind(micro_2_micro_1_5kb,some_count)
}
c42b_micro_2_5kb$Micro_C_1_Billion[c42b_micro_2_5kb$loopname %in% micro_2_micro_1_5kb$loopname] <- 1

micro_2_micro_3_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_2_micro_3_5kb <-rbind(micro_2_micro_3_5kb,some_count)
}
c42b_micro_2_5kb$Micro_C_3_Billion[c42b_micro_2_5kb$loopname %in% micro_2_micro_3_5kb$loopname] <- 1



micro_3_hic_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_5kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_3_hic_5kb <-rbind(micro_3_hic_5kb,some_count)
}
c42b_micro_3_5kb$Hi_C[c42b_micro_3_5kb$loopname %in% micro_3_hic_5kb$loopname] <- 1

micro_3_micro_1_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_3_micro_1_5kb <-rbind(micro_3_micro_1_5kb,some_count)
}
c42b_micro_3_5kb$Micro_C_1_Billion[c42b_micro_3_5kb$loopname %in% micro_3_micro_1_5kb$loopname] <- 1

micro_3_micro_2_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_3_micro_2_5kb <-rbind(micro_3_micro_2_5kb,some_count)
}
c42b_micro_3_5kb$Micro_C_2_Billion[c42b_micro_3_5kb$loopname %in% micro_3_micro_2_5kb$loopname] <- 1



hic_micro_3_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  hic_micro_3_5kb <-rbind(hic_micro_3_5kb,some_count)
}
c42b_hic_5kb$Micro_C_3_Billion[c42b_hic_5kb$loopname %in% hic_micro_3_5kb$loopname] <- 1

hic_micro_1_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  hic_micro_1_5kb <-rbind(hic_micro_1_5kb,some_count)
}
c42b_hic_5kb$Micro_C_1_Billion[c42b_hic_5kb$loopname %in% hic_micro_1_5kb$loopname] <- 1

hic_micro_2_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  hic_micro_2_5kb <-rbind(hic_micro_2_5kb,some_count)
}
c42b_hic_5kb$Micro_C_2_Billion[c42b_hic_5kb$loopname %in% hic_micro_2_5kb$loopname] <- 1





c42b_hic_1kb <-  read.table("HiC_Mustache-1kb-all-loop.tsv")
c42b_micro_1_1kb <- read.table("1_Billion_Mustache-1kb-all-loop.tsv")
c42b_micro_3_1kb <- read.table("3_billion_high_Mustache-1kb-all-loop.tsv")
c42b_micro_2_1kb <- read.table("2_Billion_Mustache-1kb-all-loop.tsv")

c42b_micro_2_1kb$loopname <- paste("c42b_2_billion",rownames(c42b_micro_2_1kb),sep = "_")
c42b_hic_1kb$loopname <- paste("c42b_2_billion",rownames(c42b_hic_1kb),sep = "_")
c42b_micro_1_1kb$loopname <- paste("c42b_1_billion",rownames(c42b_micro_1_1kb),sep = "_")
c42b_micro_3_1kb$loopname <- paste("c42b_3_billion",rownames(c42b_micro_3_1kb),sep = "_")

c42b_hic_1kb <- c42b_hic_1kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 1, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0)
c42b_micro_1_1kb <- c42b_micro_1_1kb %>% mutate(Micro_C_1_Billion = 1, Hi_C = 0, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0, Chicago = 0)
c42b_micro_3_1kb <- c42b_micro_3_1kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 0, Micro_C_2_Billion = 0,Micro_C_3_Billion = 1, Chicago = 0)
c42b_micro_2_1kb <- c42b_micro_2_1kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 0, Micro_C_2_Billion = 1,Micro_C_3_Billion = 0, Chicago = 0)

c42b_micro_2_1kb$V1 <- gsub("chr","",c42b_micro_2_1kb$V1)
c42b_micro_2_1kb$V4 <- gsub("chr","",c42b_micro_2_1kb$V4)


x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")


micro_1_hic_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_1kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  micro_1_hic_1kb <-rbind(micro_1_hic_1kb,some_count)
}
c42b_micro_1_1kb$Hi_C[c42b_micro_1_1kb$loopname %in% micro_1_hic_1kb$loopname] <- 1

micro_1_micro_2_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  micro_1_micro_2_1kb <-rbind(micro_1_micro_2_1kb,some_count)
}
c42b_micro_1_1kb$Micro_C_2_Billion[c42b_micro_1_1kb$loopname %in% micro_1_micro_2_1kb$loopname] <- 1

micro_1_micro_3_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  micro_1_micro_3_1kb <-rbind(micro_1_micro_3_1kb,some_count)
}
c42b_micro_1_1kb$Micro_C_3_Billion[c42b_micro_1_1kb$loopname %in% micro_1_micro_3_1kb$loopname] <- 1



micro_2_hic_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_1kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  micro_2_hic_1kb <-rbind(micro_2_hic_1kb,some_count)
}
c42b_micro_2_1kb$Hi_C[c42b_micro_2_1kb$loopname %in% micro_2_hic_1kb$loopname] <- 1

micro_2_micro_1_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  micro_2_micro_1_1kb <-rbind(micro_2_micro_1_1kb,some_count)
}
c42b_micro_2_1kb$Micro_C_1_Billion[c42b_micro_2_1kb$loopname %in% micro_2_micro_1_1kb$loopname] <- 1

micro_2_micro_3_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  micro_2_micro_3_1kb <-rbind(micro_2_micro_3_1kb,some_count)
}
c42b_micro_2_1kb$Micro_C_3_Billion[c42b_micro_2_1kb$loopname %in% micro_2_micro_3_1kb$loopname] <- 1



micro_3_hic_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_1kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  micro_3_hic_1kb <-rbind(micro_3_hic_1kb,some_count)
}
c42b_micro_3_1kb$Hi_C[c42b_micro_3_1kb$loopname %in% micro_3_hic_1kb$loopname] <- 1

micro_3_micro_1_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  micro_3_micro_1_1kb <-rbind(micro_3_micro_1_1kb,some_count)
}
c42b_micro_3_1kb$Micro_C_1_Billion[c42b_micro_3_1kb$loopname %in% micro_3_micro_1_1kb$loopname] <- 1

micro_3_micro_2_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  micro_3_micro_2_1kb <-rbind(micro_3_micro_2_1kb,some_count)
}
c42b_micro_3_1kb$Micro_C_2_Billion[c42b_micro_3_1kb$loopname %in% micro_3_micro_2_1kb$loopname] <- 1



hic_micro_3_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  hic_micro_3_1kb <-rbind(hic_micro_3_1kb,some_count)
}
c42b_hic_1kb$Micro_C_3_Billion[c42b_hic_1kb$loopname %in% hic_micro_3_1kb$loopname] <- 1

hic_micro_1_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  hic_micro_1_1kb <-rbind(hic_micro_1_1kb,some_count)
}
c42b_hic_1kb$Micro_C_1_Billion[c42b_hic_1kb$loopname %in% hic_micro_1_1kb$loopname] <- 1

hic_micro_2_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 2000, by = c("V2","V3","V5","V6"))
  hic_micro_2_1kb <-rbind(hic_micro_2_1kb,some_count)
}
c42b_hic_1kb$Micro_C_2_Billion[c42b_hic_1kb$loopname %in% hic_micro_2_1kb$loopname] <- 1




c42b_hic_2kb <-  read.table("HiC_Mustache-2kb-all-loop.tsv")
c42b_micro_1_2kb <- read.table("1_Billion_Mustache-2kb-all-loop.tsv")
c42b_micro_3_2kb <- read.table("3_billion_high_Mustache-2kb-all-loop.tsv")
c42b_micro_2_2kb <- read.table("2_Billion_Mustache-2kb-all-loop.tsv")

c42b_micro_2_2kb$loopname <- paste("c42b_2_billion",rownames(c42b_micro_2_2kb),sep = "_")
c42b_hic_2kb$loopname <- paste("c42b_2_billion",rownames(c42b_hic_2kb),sep = "_")
c42b_micro_1_2kb$loopname <- paste("c42b_1_billion",rownames(c42b_micro_1_2kb),sep = "_")
c42b_micro_3_2kb$loopname <- paste("c42b_3_billion",rownames(c42b_micro_3_2kb),sep = "_")

c42b_hic_2kb <- c42b_hic_2kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 1, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0)
c42b_micro_1_2kb <- c42b_micro_1_2kb %>% mutate(Micro_C_1_Billion = 1, Hi_C = 0, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0, Chicago = 0)
c42b_micro_3_2kb <- c42b_micro_3_2kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 0, Micro_C_2_Billion = 0,Micro_C_3_Billion = 1, Chicago = 0)
c42b_micro_2_2kb <- c42b_micro_2_2kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 0, Micro_C_2_Billion = 1,Micro_C_3_Billion = 0, Chicago = 0)

c42b_micro_2_2kb$V1 <- gsub("chr","",c42b_micro_2_2kb$V1)
c42b_micro_2_2kb$V4 <- gsub("chr","",c42b_micro_2_2kb$V4)

x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")


micro_1_hic_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_2kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_1_hic_2kb <-rbind(micro_1_hic_2kb,some_count)
}
c42b_micro_1_2kb$Hi_C[c42b_micro_1_2kb$loopname %in% micro_1_hic_2kb$loopname] <- 1

micro_1_micro_2_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_1_micro_2_2kb <-rbind(micro_1_micro_2_2kb,some_count)
}
c42b_micro_1_2kb$Micro_C_2_Billion[c42b_micro_1_2kb$loopname %in% micro_1_micro_2_2kb$loopname] <- 1

micro_1_micro_3_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_1_micro_3_2kb <-rbind(micro_1_micro_3_2kb,some_count)
}
c42b_micro_1_2kb$Micro_C_3_Billion[c42b_micro_1_2kb$loopname %in% micro_1_micro_3_2kb$loopname] <- 1



micro_2_hic_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_2kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_2_hic_2kb <-rbind(micro_2_hic_2kb,some_count)
}
c42b_micro_2_2kb$Hi_C[c42b_micro_2_2kb$loopname %in% micro_2_hic_2kb$loopname] <- 1

micro_2_micro_1_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_2_micro_1_2kb <-rbind(micro_2_micro_1_2kb,some_count)
}
c42b_micro_2_2kb$Micro_C_1_Billion[c42b_micro_2_2kb$loopname %in% micro_2_micro_1_2kb$loopname] <- 1

micro_2_micro_3_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_2_micro_3_2kb <-rbind(micro_2_micro_3_2kb,some_count)
}
c42b_micro_2_2kb$Micro_C_3_Billion[c42b_micro_2_2kb$loopname %in% micro_2_micro_3_2kb$loopname] <- 1



micro_3_hic_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_2kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_3_hic_2kb <-rbind(micro_3_hic_2kb,some_count)
}
c42b_micro_3_2kb$Hi_C[c42b_micro_3_2kb$loopname %in% micro_3_hic_2kb$loopname] <- 1

micro_3_micro_1_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_3_micro_1_2kb <-rbind(micro_3_micro_1_2kb,some_count)
}
c42b_micro_3_2kb$Micro_C_1_Billion[c42b_micro_3_2kb$loopname %in% micro_3_micro_1_2kb$loopname] <- 1

micro_3_micro_2_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_3_micro_2_2kb <-rbind(micro_3_micro_2_2kb,some_count)
}
c42b_micro_3_2kb$Micro_C_2_Billion[c42b_micro_3_2kb$loopname %in% micro_3_micro_2_2kb$loopname] <- 1



hic_micro_3_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  hic_micro_3_2kb <-rbind(hic_micro_3_2kb,some_count)
}
c42b_hic_2kb$Micro_C_3_Billion[c42b_hic_2kb$loopname %in% hic_micro_3_2kb$loopname] <- 1

hic_micro_1_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  hic_micro_1_2kb <-rbind(hic_micro_1_2kb,some_count)
}
c42b_hic_2kb$Micro_C_1_Billion[c42b_hic_2kb$loopname %in% hic_micro_1_2kb$loopname] <- 1

hic_micro_2_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  hic_micro_2_2kb <-rbind(hic_micro_2_2kb,some_count)
}
c42b_hic_2kb$Micro_C_2_Billion[c42b_hic_2kb$loopname %in% hic_micro_2_2kb$loopname] <- 1




c42b_hic_10kb <-  read.table("HiC_Mustache-10kb-all-loop.tsv")
c42b_micro_1_10kb <- read.table("1_Billion_Mustache-10kb-all-loop.tsv")
c42b_micro_3_10kb <- read.table("3_billion_low_Mustache-10kb-all-loop.tsv")
c42b_micro_2_10kb <- read.table("2_Billion_Mustache-10kb-all-loop.tsv")

c42b_micro_2_10kb$loopname <- paste("c42b_2_billion",rownames(c42b_micro_2_10kb),sep = "_")
c42b_hic_10kb$loopname <- paste("c42b_2_billion",rownames(c42b_hic_10kb),sep = "_")
c42b_micro_1_10kb$loopname <- paste("c42b_1_billion",rownames(c42b_micro_1_10kb),sep = "_")
c42b_micro_3_10kb$loopname <- paste("c42b_3_billion",rownames(c42b_micro_3_10kb),sep = "_")

c42b_hic_10kb <- c42b_hic_10kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 1, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0)
c42b_micro_1_10kb <- c42b_micro_1_10kb %>% mutate(Micro_C_1_Billion = 1, Hi_C = 0, Micro_C_2_Billion = 0,Micro_C_3_Billion = 0, Chicago = 0)
c42b_micro_3_10kb <- c42b_micro_3_10kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 0, Micro_C_2_Billion = 0,Micro_C_3_Billion = 1, Chicago = 0)
c42b_micro_2_10kb <- c42b_micro_2_10kb %>% mutate(Micro_C_1_Billion = 0, Hi_C = 0, Micro_C_2_Billion = 1,Micro_C_3_Billion = 0, Chicago = 0)

c42b_micro_2_10kb$V1 <- gsub("chr","",c42b_micro_2_10kb$V1)
c42b_micro_2_10kb$V4 <- gsub("chr","",c42b_micro_2_10kb$V4)

x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")


micro_1_hic_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_10kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_1_hic_10kb <-rbind(micro_1_hic_10kb,some_count)
}
c42b_micro_1_10kb$Hi_C[c42b_micro_1_10kb$loopname %in% micro_1_hic_10kb$loopname] <- 1

micro_1_micro_2_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_1_micro_2_10kb <-rbind(micro_1_micro_2_10kb,some_count)
}
c42b_micro_1_10kb$Micro_C_2_Billion[c42b_micro_1_10kb$loopname %in% micro_1_micro_2_10kb$loopname] <- 1

micro_1_micro_3_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_1_micro_3_10kb <-rbind(micro_1_micro_3_10kb,some_count)
}
c42b_micro_1_10kb$Micro_C_3_Billion[c42b_micro_1_10kb$loopname %in% micro_1_micro_3_10kb$loopname] <- 1



micro_2_hic_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_10kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_2_hic_10kb <-rbind(micro_2_hic_10kb,some_count)
}
c42b_micro_2_10kb$Hi_C[c42b_micro_2_10kb$loopname %in% micro_2_hic_10kb$loopname] <- 1

micro_2_micro_1_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_2_micro_1_10kb <-rbind(micro_2_micro_1_10kb,some_count)
}
c42b_micro_2_10kb$Micro_C_1_Billion[c42b_micro_2_10kb$loopname %in% micro_2_micro_1_10kb$loopname] <- 1

micro_2_micro_3_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_2_micro_3_10kb <-rbind(micro_2_micro_3_10kb,some_count)
}
c42b_micro_2_10kb$Micro_C_3_Billion[c42b_micro_2_10kb$loopname %in% micro_2_micro_3_10kb$loopname] <- 1



micro_3_hic_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_10kb %>% filter(V1 == i)
  chr_reference <- c42b_hic_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_3_hic_10kb <-rbind(micro_3_hic_10kb,some_count)
}
c42b_micro_3_10kb$Hi_C[c42b_micro_3_10kb$loopname %in% micro_3_hic_10kb$loopname] <- 1

micro_3_micro_1_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_3_micro_1_10kb <-rbind(micro_3_micro_1_10kb,some_count)
}
c42b_micro_3_10kb$Micro_C_1_Billion[c42b_micro_3_10kb$loopname %in% micro_3_micro_1_10kb$loopname] <- 1

micro_3_micro_2_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_3_micro_2_10kb <-rbind(micro_3_micro_2_10kb,some_count)
}
c42b_micro_3_10kb$Micro_C_2_Billion[c42b_micro_3_10kb$loopname %in% micro_3_micro_2_10kb$loopname] <- 1



hic_micro_3_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  hic_micro_3_10kb <-rbind(hic_micro_3_10kb,some_count)
}
c42b_hic_10kb$Micro_C_3_Billion[c42b_hic_10kb$loopname %in% hic_micro_3_10kb$loopname] <- 1

hic_micro_1_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  hic_micro_1_10kb <-rbind(hic_micro_1_10kb,some_count)
}
c42b_hic_10kb$Micro_C_1_Billion[c42b_hic_10kb$loopname %in% hic_micro_1_10kb$loopname] <- 1

hic_micro_2_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_hic_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  hic_micro_2_10kb <-rbind(hic_micro_2_10kb,some_count)
}
c42b_hic_10kb$Micro_C_2_Billion[c42b_hic_10kb$loopname %in% hic_micro_2_10kb$loopname] <- 1



micro_1_chicago_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_5kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_1_chicago_5kb <-rbind(micro_1_chicago_5kb,some_count)
}
c42b_micro_1_5kb$Chicago[c42b_micro_1_5kb$loopname %in% micro_1_chicago_5kb$loopname] <- 1

micro_2_chicago_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_5kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_2_chicago_5kb <-rbind(micro_2_chicago_5kb,some_count)
}
c42b_micro_2_5kb$Chicago[c42b_micro_2_5kb$loopname %in% micro_2_chicago_5kb$loopname] <- 1

micro_3_chicago_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_5kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_3_chicago_5kb <-rbind(micro_3_chicago_5kb,some_count)
}
c42b_micro_3_5kb$Chicago[c42b_micro_3_5kb$loopname %in% micro_3_chicago_5kb$loopname] <- 1


chicago_micro_3_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  chicago_micro_3_5kb <-rbind(chicago_micro_3_5kb,some_count)
}
c42b_chicago_5kb$Micro_C_3_Billion[c42b_chicago_5kb$loopname %in% chicago_micro_3_5kb$loopname] <- 1

chicago_micro_1_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  chicago_micro_1_5kb <-rbind(chicago_micro_1_5kb,some_count)
}
c42b_chicago_5kb$Micro_C_1_Billion[c42b_chicago_5kb$loopname %in% chicago_micro_1_5kb$loopname] <- 1

chicago_micro_2_5kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_5kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_5kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  chicago_micro_2_5kb <-rbind(chicago_micro_2_5kb,some_count)
}
c42b_chicago_5kb$Micro_C_2_Billion[c42b_chicago_5kb$loopname %in% chicago_micro_2_5kb$loopname] <- 1


micro_1_chicago_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_1kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_1_chicago_1kb <-rbind(micro_1_chicago_1kb,some_count)
}
c42b_micro_1_1kb$Chicago[c42b_micro_1_1kb$loopname %in% micro_1_chicago_1kb$loopname] <- 1

micro_2_chicago_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_1kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_2_chicago_1kb <-rbind(micro_2_chicago_1kb,some_count)
}
c42b_micro_2_1kb$Chicago[c42b_micro_2_1kb$loopname %in% micro_2_chicago_1kb$loopname] <- 1

micro_3_chicago_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_1kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  micro_3_chicago_1kb <-rbind(micro_3_chicago_1kb,some_count)
}
c42b_micro_3_1kb$Chicago[c42b_micro_3_1kb$loopname %in% micro_3_chicago_1kb$loopname] <- 1

chicago_micro_3_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  chicago_micro_3_1kb <-rbind(chicago_micro_3_1kb,some_count)
}
c42b_chicago_1kb$Micro_C_3_Billion[c42b_chicago_1kb$loopname %in% chicago_micro_3_1kb$loopname] <- 1

chicago_micro_1_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  chicago_micro_1_1kb <-rbind(chicago_micro_1_1kb,some_count)
}
c42b_chicago_1kb$Micro_C_1_Billion[c42b_chicago_1kb$loopname %in% chicago_micro_1_1kb$loopname] <- 1

chicago_micro_2_1kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_1kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_1kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 10000, by = c("V2","V3","V5","V6"))
  chicago_micro_2_1kb <-rbind(chicago_micro_2_1kb,some_count)
}
c42b_chicago_1kb$Micro_C_2_Billion[c42b_chicago_1kb$loopname %in% chicago_micro_2_1kb$loopname] <- 1

micro_1_chicago_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_2kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_1_chicago_2kb <-rbind(micro_1_chicago_2kb,some_count)
}
c42b_micro_1_2kb$Chicago[c42b_micro_1_2kb$loopname %in% micro_1_chicago_2kb$loopname] <- 1

micro_2_chicago_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_2kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_2_chicago_2kb <-rbind(micro_2_chicago_2kb,some_count)
}
c42b_micro_2_2kb$Chicago[c42b_micro_2_2kb$loopname %in% micro_2_chicago_2kb$loopname] <- 1

micro_3_chicago_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_2kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  micro_3_chicago_2kb <-rbind(micro_3_chicago_2kb,some_count)
}
c42b_micro_3_2kb$Chicago[c42b_micro_3_2kb$loopname %in% micro_3_chicago_2kb$loopname] <- 1

chicago_micro_3_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  chicago_micro_3_2kb <-rbind(chicago_micro_3_2kb,some_count)
}
c42b_chicago_2kb$Micro_C_3_Billion[c42b_chicago_2kb$loopname %in% chicago_micro_3_2kb$loopname] <- 1

chicago_micro_1_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  chicago_micro_1_2kb <-rbind(chicago_micro_1_2kb,some_count)
}
c42b_chicago_2kb$Micro_C_1_Billion[c42b_chicago_2kb$loopname %in% chicago_micro_1_2kb$loopname] <- 1

chicago_micro_2_2kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_2kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_2kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 4000, by = c("V2","V3","V5","V6"))
  chicago_micro_2_2kb <-rbind(chicago_micro_2_2kb,some_count)
}
c42b_chicago_2kb$Micro_C_2_Billion[c42b_chicago_2kb$loopname %in% chicago_micro_2_2kb$loopname] <- 1


micro_1_chicago_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_1_10kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_1_chicago_10kb <-rbind(micro_1_chicago_10kb,some_count)
}
c42b_micro_1_10kb$Chicago[c42b_micro_1_10kb$loopname %in% micro_1_chicago_10kb$loopname] <- 1

micro_2_chicago_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_2_10kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_2_chicago_10kb <-rbind(micro_2_chicago_10kb,some_count)
}
c42b_micro_2_10kb$Chicago[c42b_micro_2_10kb$loopname %in% micro_2_chicago_10kb$loopname] <- 1

micro_3_chicago_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_micro_3_10kb %>% filter(V1 == i)
  chr_reference <- c42b_chicago_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  micro_3_chicago_10kb <-rbind(micro_3_chicago_10kb,some_count)
}
c42b_micro_3_10kb$Chicago[c42b_micro_3_10kb$loopname %in% micro_3_chicago_10kb$loopname] <- 1

chicago_micro_3_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_3_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  chicago_micro_3_10kb <-rbind(chicago_micro_3_10kb,some_count)
}
c42b_chicago_10kb$Micro_C_3_Billion[c42b_chicago_10kb$loopname %in% chicago_micro_3_10kb$loopname] <- 1

chicago_micro_1_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_1_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  chicago_micro_1_10kb <-rbind(chicago_micro_1_10kb,some_count)
}
c42b_chicago_10kb$Micro_C_1_Billion[c42b_chicago_10kb$loopname %in% chicago_micro_1_10kb$loopname] <- 1

chicago_micro_2_10kb <- data.frame()
for (i in x){
  chr_template <- c42b_chicago_10kb %>% filter(V1 == i)
  chr_reference <- c42b_micro_2_10kb %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= 20000, by = c("V2","V3","V5","V6"))
  chicago_micro_2_10kb <-rbind(chicago_micro_2_10kb,some_count)
}
c42b_chicago_10kb$Micro_C_2_Billion[c42b_chicago_10kb$loopname %in% chicago_micro_2_10kb$loopname] <- 1







#Add the chr back to the loops for future analysis in case needed.
c42b_hic_5kb$V1 <- paste0("chr",c42b_hic_5kb$V1)
c42b_hic_5kb$V4 <- paste0("chr",c42b_hic_5kb$V4)

c42b_hic_1kb$V1 <- paste0("chr",c42b_hic_1kb$V1)
c42b_hic_1kb$V4 <- paste0("chr",c42b_hic_1kb$V4)

c42b_hic_2kb$V1 <- paste0("chr",c42b_hic_2kb$V1)
c42b_hic_2kb$V4 <- paste0("chr",c42b_hic_2kb$V4)

c42b_hic_10kb$V1 <- paste0("chr",c42b_hic_10kb$V1)
c42b_hic_10kb$V4 <- paste0("chr",c42b_hic_10kb$V4)

c42b_micro_1_5kb$V1 <- paste0("chr",c42b_micro_1_5kb$V1)
c42b_micro_1_5kb$V4 <- paste0("chr",c42b_micro_1_5kb$V4)

c42b_micro_1_1kb$V1 <- paste0("chr",c42b_micro_1_1kb$V1)
c42b_micro_1_1kb$V4 <- paste0("chr",c42b_micro_1_1kb$V4)

c42b_micro_1_2kb$V1 <- paste0("chr",c42b_micro_1_2kb$V1)
c42b_micro_1_2kb$V4 <- paste0("chr",c42b_micro_1_2kb$V4)

c42b_micro_1_10kb$V1 <- paste0("chr",c42b_micro_1_10kb$V1)
c42b_micro_1_10kb$V4 <- paste0("chr",c42b_micro_1_10kb$V4)

c42b_micro_2_5kb$V1 <- paste0("chr",c42b_micro_2_5kb$V1)
c42b_micro_2_5kb$V4 <- paste0("chr",c42b_micro_2_5kb$V4)

c42b_micro_2_1kb$V1 <- paste0("chr",c42b_micro_2_1kb$V1)
c42b_micro_2_1kb$V4 <- paste0("chr",c42b_micro_2_1kb$V4)

c42b_micro_2_2kb$V1 <- paste0("chr",c42b_micro_2_2kb$V1)
c42b_micro_2_2kb$V4 <- paste0("chr",c42b_micro_2_2kb$V4)

c42b_micro_2_10kb$V1 <- paste0("chr",c42b_micro_2_10kb$V1)
c42b_micro_2_10kb$V4 <- paste0("chr",c42b_micro_2_10kb$V4)

c42b_micro_3_5kb$V1 <- paste0("chr",c42b_micro_3_5kb$V1)
c42b_micro_3_5kb$V4 <- paste0("chr",c42b_micro_3_5kb$V4)

c42b_micro_3_1kb$V1 <- paste0("chr",c42b_micro_3_1kb$V1)
c42b_micro_3_1kb$V4 <- paste0("chr",c42b_micro_3_1kb$V4)

c42b_micro_3_2kb$V1 <- paste0("chr",c42b_micro_3_2kb$V1)
c42b_micro_3_2kb$V4 <- paste0("chr",c42b_micro_3_2kb$V4)

c42b_micro_3_10kb$V1 <- paste0("chr",c42b_micro_3_10kb$V1)
c42b_micro_3_10kb$V4 <- paste0("chr",c42b_micro_3_10kb$V4)

c42b_chicago_5kb$V1 <- paste0("chr",c42b_chicago_5kb$V1)
c42b_chicago_5kb$V4 <- paste0("chr",c42b_chicago_5kb$V4)

c42b_chicago_1kb$V1 <- paste0("chr",c42b_chicago_1kb$V1)
c42b_chicago_1kb$V4 <- paste0("chr",c42b_chicago_1kb$V4)

c42b_chicago_2kb$V1 <- paste0("chr",c42b_chicago_2kb$V1)
c42b_chicago_2kb$V4 <- paste0("chr",c42b_chicago_2kb$V4)

c42b_chicago_10kb$V1 <- paste0("chr",c42b_chicago_10kb$V1)
c42b_chicago_10kb$V4 <- paste0("chr",c42b_chicago_10kb$V4)


#Select necessary informations only, then export them as a table.
c42b_hic_score_5kb <- c42b_hic_5kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion)
c42b_hic_score_1kb <- c42b_hic_1kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion)
c42b_hic_score_2kb <- c42b_hic_2kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion)
c42b_hic_score_10kb <- c42b_hic_10kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion)

c42b_microc1_score_10kb <- c42b_micro_1_10kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_microc1_score_1kb <- c42b_micro_1_1kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_microc1_score_2kb <- c42b_micro_1_2kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_microc1_score_5kb <- c42b_micro_1_5kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)

c42b_microc2_score_10kb <- c42b_micro_2_10kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_microc2_score_1kb <- c42b_micro_2_1kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_microc2_score_2kb <- c42b_micro_2_2kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_microc2_score_5kb <- c42b_micro_2_5kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)

c42b_microc3_score_10kb <- c42b_micro_3_10kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_microc3_score_1kb <- c42b_micro_3_1kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_microc3_score_2kb <- c42b_micro_3_2kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_microc3_score_5kb <- c42b_micro_3_5kb %>% select(V1,V2,V3,V4,V5,V6,Hi_C ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)

c42b_chicago_score_10kb <- c42b_chicago_10kb %>% select(V1,V2,V3,V4,V5,V6,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_chicago_score_1kb <- c42b_chicago_1kb %>% select(V1,V2,V3,V4,V5,V6 ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_chicago_score_2kb <- c42b_chicago_2kb %>% select(V1,V2,V3,V4,V5,V6,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)
c42b_chicago_score_5kb <- c42b_chicago_5kb %>% select(V1,V2,V3,V4,V5,V6 ,Micro_C_1_Billion, Micro_C_2_Billion, Micro_C_3_Billion,Chicago)


write.table(c42b_hic_score_5kb, "HiC_loop_sharing_with_score_5kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_hic_score_1kb, "HiC_loop_sharing_with_score_1kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_hic_score_2kb, "HiC_loop_sharing_with_score_2kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_hic_score_10kb, "HiC_loop_sharing_with_score_10kb.csv", row.names = FALSE, quote = FALSE, sep ="," )

write.table(c42b_microc1_score_5kb, "Micro_C_1_Billion_loop_sharing_with_score_5kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_microc1_score_1kb, "Micro_C_1_Billion_loop_sharing_with_score_1kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_microc1_score_2kb, "Micro_C_1_Billion_loop_sharing_with_score_2kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_microc1_score_10kb, "Micro_C_1_Billion_loop_sharing_with_score_10kb.csv", row.names = FALSE, quote = FALSE, sep ="," )

write.table(c42b_microc2_score_5kb, "Micro_C_2_Billion_loop_sharing_with_score_5kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_microc2_score_1kb, "Micro_C_2_Billion_loop_sharing_with_score_1kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_microc2_score_2kb, "Micro_C_2_Billion_loop_sharing_with_score_2kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_microc2_score_10kb, "Micro_C_2_Billion_loop_sharing_with_score_10kb.csv", row.names = FALSE, quote = FALSE, sep ="," )

write.table(c42b_microc3_score_5kb, "Micro_C_3_Billion_loop_sharing_with_score_5kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_microc3_score_1kb, "Micro_C_3_Billion_loop_sharing_with_score_1kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_microc3_score_2kb, "Micro_C_3_Billion_loop_sharing_with_score_2kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_microc3_score_10kb, "Micro_C_3_Billion_loop_sharing_with_score_10kb.csv", row.names = FALSE, quote = FALSE, sep ="," )

write.table(c42b_chicago_score_5kb, "Chicago_loop_sharing_with_score_5kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_chicago_score_1kb, "Chicago_loop_sharing_with_score_1kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_chicago_score_2kb, "Chicago_loop_sharing_with_score_2kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
write.table(c42b_chicago_score_10kb, "Chicago_loop_sharing_with_score_10kb.csv", row.names = FALSE, quote = FALSE, sep ="," )
