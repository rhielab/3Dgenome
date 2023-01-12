library(dplyr)
library(fuzzyjoin)

topdom_comparison <- function(
topdom_bed1,
topdom_bed2,
label_1 = "topdom_bed1",
label_2 = "topdom_bed2",
resolution
){
topdom_bed1$V1 <- gsub("chr","",topdom_bed1$V1)
topdom_bed2$V1 <- gsub("chr","",topdom_bed2$V1)

#put loop name
topdom_bed1$loopname <- paste(label_1,rownames(topdom_bed1),sep="_")
topdom_bed2$loopname <- paste(label_2,rownames(topdom_bed2),sep="_")

#put comparison column
topdom_bed1$comparison <- 0
topdom_bed2$comparison <- 0


#only perform analysis in "domain" within V4 column
topdom_bed1 <- topdom_bed1 %>% filter(V4 == "domain")
topdom_bed2 <- topdom_bed2 %>% filter(V4 == "domain")

x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")

bed1_to_bed2 <- data.frame()
bed2_to_bed1 <- data.frame()
for (i in x){
  chr_template <- topdom_bed1 %>% filter(V1 == i)
  chr_reference <- topdom_bed2 %>% filter(V1 == i)
  some_count <- difference_semi_join(chr_template,chr_reference, max_dist= resolution*2, by = c("V2","V3"))
  some_count2 <- difference_semi_join(chr_reference,chr_template, max_dist= resolution*2, by = c("V2","V3"))
  bed1_to_bed2 <- rbind(bed1_to_bed2,some_count)
  bed2_to_bed1 <- rbind(bed2_to_bed1,some_count2)
}

#update comparison matrix
bed1_to_bed2$comparison <- 1
bed2_to_bed1$comparison <- 1
topdom_bed1$comparison[topdom_bed1$loopname %in% bed1_to_bed2$loopname] <- 1
topdom_bed2$comparison[topdom_bed2$loopname %in% bed2_to_bed1$loopname] <- 1

#change column comparison to the labels
colnames(topdom_bed1)[ncol(topdom_bed1)] <- label_2
colnames(bed1_to_bed2)[ncol(bed1_to_bed2)] <- label_2
colnames(topdom_bed2)[ncol(topdom_bed2)] <- label_1
colnames(bed2_to_bed1)[ncol(bed2_to_bed1)] <- label_1

#create list to save wanted outputs
final <- list()
final[[1]] <- bed1_to_bed2
final[[2]] <- bed2_to_bed1
final[[3]] <- topdom_bed1
final[[4]] <- topdom_bed2

#remove loopname column, change column comparison to the labels, put chr back in V1
final <- lapply(final, function(x){
	x <- subset(x, select = -loopname)
	x$V1 <- paste0('chr', x$V1)
	return(x)
 	}) #end lapply

return(final)
} #end function
