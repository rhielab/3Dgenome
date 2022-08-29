library(dplyr)
library(fuzzyjoin)

check_REGelement <- function(
loop_data,
reg_data,
reg_type,
output_path,
resolution,
prefix
){
#modify the data
names(loop_data)[1:6] <- c('V1','V2','V3','V4','V5','V6')
reg_data$name <- paste0(reg_type, rownames(reg_data))
names(reg_data)[1:3] <- c('V1','V2','V3')

#Check overlap with reg
loop_data$loopname <- paste(reg_type,rownames(loop_data),sep = "_")

loop_data$V1 <- gsub("chr","",loop_data$V1)
loop_data$V4 <- gsub("chr","",loop_data$V4)
reg_data$V1 <- gsub("chr","",reg_data$V1)

x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")

anchor1 <- loop_data %>% select(V1,V2,V3,loopname)
anchor2 <- loop_data %>% select(V4,V5,V6,loopname)

colnames(anchor2)[1:3] <- c("V1","V2","V3")

#For regulatory element
anchor_reg_data <- data.frame()
for (i in x){
        chr_fithic <- reg_data %>% filter(V1 == i)
        #For anchor1
  chr_mustache <- anchor1 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic, chr_mustache, max_dist= 2*resolution, by = c("V2","V3"))
  anchor_reg_data <- rbind(anchor_reg_data,some_count)

  #For anchor2
  chr_mustache <- anchor2 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_fithic, chr_mustache, max_dist= 2*resolution, by = c("V2","V3"))
  anchor_reg_data <- rbind(anchor_reg_data,some_count)
        }

#Give chr back to V1
anchor_reg_data$V1 <- paste0('chr', anchor_reg_data$V1)
reg_data$V1 <- paste0('chr', reg_data$V1)

reg_inLoop <- unique(anchor_reg_data)
reg_nonLoop <- reg_data[!(reg_data$name %in% reg_inLoop$name),]

summary <- t(data.frame(
                N_reg = nrow(reg_data),
                N_reg_inLoop = nrow(reg_inLoop),
                N_reg_nonLoop = nrow(reg_nonLoop),
                PCT_reg_inLoop = nrow(reg_inLoop)/nrow(reg_data) ,
                PCT_reg_nonLoop = nrow(reg_nonLoop)/nrow(reg_data)
                ))
rownames(summary) <- gsub('reg',reg_type, rownames(summary))
#output files
write.table(reg_inLoop, file = paste0(output_path, '/', prefix,'_inLoop.bed'), quote = F, row.names = F, col.names = F, sep = '\t')
write.table(reg_nonLoop, file = paste0(output_path, '/', prefix,'_nonLoop.bed'), quote = F, row.names = F, col.names = F, sep = '\t')
write.table(summary, file = paste0(output_path, '/', prefix,'_summary.txt'), quote = F, row.names = T, col.names = F, sep = '\t')
}
