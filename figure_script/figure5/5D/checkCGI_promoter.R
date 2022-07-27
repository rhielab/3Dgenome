library(dplyr)
library(fuzzyjoin)

#Check CIP content for in-loop promoter and non in-loop promoter

checkCGI_promoter <- function(
loop_data,
promoter_data_path = '/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/figure4_loop_reg/input/C42B_TSS_FPKM_greater_than_0.5.tsv',
CGI_reference_path = '/project/rhie_130/suhn/zexunwu/Micro_C_paper_BE/script_testing/check_CGIpromoter/input/REF/hg38.gencode.CGI_UCSCbrowser_full.txt',
output_path,
prefix
){
#Load the data
promoter_data <- read.table(file = promoter_data_path)
CGI_reference <- read.table(file = CGI_reference_path)
promoter_data$name <- paste0('promoter_', rownames(promoter_data))
names(CGI_reference)[1:3] <- c('V1','V2','V3')
names(promoter_data)[1:3] <- c('V1','V2','V3')

#Check overlap with promoter
loop_data$loopname <- paste(prefix,rownames(loop_data),sep = "_")

loop_data$V1 <- gsub("chr","",loop_data$V1)
loop_data$V4 <- gsub("chr","",loop_data$V4)
promoter_data$V1 <- gsub("chr","",promoter_data$V1)

x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")

anchor1 <- loop_data %>% select(V1,V2,V3,loopname)
anchor2 <- loop_data %>% select(V4,V5,V6,loopname)

colnames(anchor2)[1] <- "V1"
colnames(anchor2)[2] <- "V2"
colnames(anchor2)[3] <- "V3"

#For promoter
anchor_promoter_data <- data.frame()
for (i in x){
        chr_fithic <- promoter_data %>% filter(V1 == i)
        #For anchor1
  chr_mustache <- anchor1 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic, chr_mustache, max_dist= 5000, by = c("V2","V3"))
  anchor_promoter_data <-rbind(anchor_promoter_data,some_count)

  #For anchor2
  chr_mustache <- anchor2 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic, chr_mustache, max_dist= 5000, by = c("V2","V3"))
  anchor_promoter_data <-rbind(anchor_promoter_data,some_count)
        }

#Give chr back to V1
anchor_promoter_data$V1 <- paste0('chr', anchor_promoter_data$V1)
promoter_data$V1 <- paste0('chr', promoter_data$V1)

promoter_inLoop <- unique(anchor_promoter_data)
promoter_nonLoop <- promoter_data[!(promoter_data$name %in% promoter_inLoop$name),]
#promoter_nonLoop <- genome_anti_join(promoter_data, promoter_inLoop, by = c('V1','V2','V3'))

CGI_promoter_inLoop <- genome_semi_join(promoter_inLoop, CGI_reference, by = c('V1','V2','V3'))
CGI_promoter_nonLoop <- genome_semi_join(promoter_nonLoop, CGI_reference, by = c('V1','V2','V3'))

summary <- t(data.frame(
                N_promoter_inLoop = nrow(promoter_inLoop),
                N_promoter_nonLoop = nrow(promoter_nonLoop),
                N_CGIpromoter_inLoop = nrow(CGI_promoter_inLoop),
                N_CGIpromoter_nonLoop =nrow(CGI_promoter_nonLoop),
                PCT_CGIpromoter_inLoop = nrow(CGI_promoter_inLoop)/nrow(promoter_inLoop) ,
                PCT_CGIpromoter_nonLoop = nrow(CGI_promoter_nonLoop)/nrow(promoter_nonLoop)
                ))
#output files
write.table(promoter_inLoop, file = paste0(output_path, '/', prefix,'_Promoter_inLoop.bed'), quote = F, row.names = F, col.names = F, sep = '\t')
write.table(promoter_nonLoop, file = paste0(output_path, '/', prefix,'_Promoter_nonLoop.bed'), quote = F, row.names = F, col.names = F, sep = '\t')
write.table(CGI_promoter_inLoop, file = paste0(output_path, '/', prefix,'_CGIPromoter_inLoop.bed'), quote = F, row.names = F, col.names = F, sep = '\t')
write.table(CGI_promoter_nonLoop, file = paste0(output_path, '/', prefix,'_CGIPromoter_nonLoop.bed'), quote = F, row.names = F, col.names = F, sep = '\t')
write.table(summary, file = paste0(output_path, '/', prefix,'_CGIPromoter_summary.txt'), quote = F, row.names = T, col.names = F, sep = '\t')
}
