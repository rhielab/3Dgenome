library(fuzzyjoin)
library(dplyr)

loop_comparison <- function(
loop_list, #Should be a list containing all the loop data with corresponding names and ordered by priority. Make sure that it get named!
priority = seq(length(loop_list)),
output_path
){
#Get the total length and names of loop list
length <- length(loop_list)
names_loop <- names(loop_list)[priority]

#reorder the loop_list based on priority
loop_list <- lapply(priority, function(x){
        return(loop_list[[x]])
        })
names(loop_list) <- names_loop

#data wrangling for all loop data in the list
loop_list <- lapply(seq(length), function(i){
        loop_data <- loop_list[[i]]
        loop_data$loopname <- paste(names_loop[i],
                                rownames(loop_data),sep = "_")

        #Create a data frame (0/1 table) for loop comparison
        loop_comp <- data.frame(matrix(0, nrow = nrow(loop_data), ncol = length))
        names(loop_comp) <- names_loop
        loop_comp[,match(names_loop[i], names(loop_comp) )] <- 1

        #merge it with loop data
        loop_data <- cbind(loop_data, loop_comp)

        #change chr column (V1 and V4)
        loop_data$V1 <- gsub("chr","",loop_data$V1)
        loop_data$V4 <- gsub("chr","",loop_data$V4)

        #done
        return(loop_data)
        }) #end lapply
names(loop_list) <- names_loop

#generate comparison table
names_comb <- t(combn(seq(length(loop_list)), m = 2))

#comparison among loop data
x <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")
for(i in seq(nrow(names_comb))){
        fir <- names_comb[i,1] #numeric
        sec <- names_comb[i,2] #numeric
        name_fir <- names_loop[fir]
        name_sec <- names_loop[sec]

        data1 <- data.frame()
        data2 <- data.frame()
        for(k in x){
                chr_1st <- loop_list[[fir]] %>% filter(V1 == k)
                chr_2nd <- loop_list[[sec]] %>% filter(V1 == k)
                some_count_1st <- difference_semi_join(chr_1st,chr_2nd, max_dist= 10000, by = c("V2","V3","V5","V6"))
                some_count_2nd <- difference_semi_join(chr_2nd,chr_1st, max_dist= 10000, by = c("V2","V3","V5","V6"))
                data1 <- rbind(data1, some_count_1st)
                data2 <- rbind(data2, some_count_2nd)
                } # end second for

        #give 1 to the category that is in comparison
        loop_list[[fir]][,match(
                                name_sec,
                                colnames(loop_list[[fir]])
                                )][ loop_list[[fir]]$loopname %in% data1$loopname ] <- 1
        loop_list[[sec]][,match(
                                name_fir,
                                colnames(loop_list[[sec]])
                                )][ loop_list[[sec]]$loopname %in% data2$loopname ] <- 1
        } #end first for

#put chr back
loop_list <- lapply(loop_list, function(x){
        x$V1 <- paste0('chr', x$V1)
        x$V4 <- paste0('chr', x$V4)
        return(x)
        }) #end lapply

#make summary table
shared_summary <- data.frame(matrix(0, nrow = length, ncol = length), row.names = names_loop)
colnames(shared_summary) <- names_loop

for(i in 1:length){
        #output the loop data with 0/1 comparison
        write.table(loop_list[[i]], file = paste0(output_path,'/',names_loop[i],'-WithShareMatrix.tsv'),
                                row.names = F, quote = F, col.names = T, sep = '\t')

        #output the share/unique loop data for each category
        summary <- numeric( length = (length - 1) )
        #step1 - make loop sharing IDs
        N_col_shareMatrix <- c( (ncol( loop_list[[i]] )-length+1):ncol( loop_list[[i]] ))
        N_col_shareMatrix <- N_col_shareMatrix[-i]
        name_loop_4comp <- names_loop[-i]
        share_ID <- apply(loop_list[[i]], 1, function(x){
                                Reduce(paste0, x[N_col_shareMatrix])
                                                }) #end apply

        #step2 - create subfolder for shared/unique data
        output_sub <- paste0(output_path,'/', names_loop[i])
        dir.create(output_sub)

        #step3 - extract shared data
        for(k in 1:(length)){
                if(k != length){
                        extract_shared_ID <- paste0('^',Reduce(paste0, rep(0, k-1)),1)
                        shared_data <- loop_list[[i]][grep(extract_shared_ID, share_ID),]
                        end_fix <- paste0('_sharingWith_', name_loop_4comp[k],'.tsv')
                        summary[k] <- nrow(shared_data)
                        } else {
                                extract_shared_ID <- paste0('^',Reduce(paste0, rep(0, k-1)))
                                shared_data <- loop_list[[i]][grep(extract_shared_ID, share_ID),]
                                end_fix <- '_unique.tsv'
                                summary <- append(summary, nrow(shared_data), i-1)
                                        }
                write.table(shared_data, file = paste0(output_sub,'/', names_loop[i], end_fix),
                                        row.names = T, quote = F, col.names = T, sep = '\t')
                }
        shared_summary[i,] <- summary
        } #end for
shared_summary$Total <- sapply(loop_list, nrow)

filecon <- file(paste0(output_path,'/loop_comparison_summary.txt'), 'w')
writeLines("#Summary table for the loop comparison analysis - Check with ROW-COLUMN combination.", con = filecon)
writeLines("#The number for each combination represents the loops in ROW shared with COLUMN.", con = filecon)
writeLines("#If ROW and COLUMN are same. The number means unique loops for ROW.", con = filecon)
write.table(shared_summary, file = filecon,
        row.names = T, quote = F, col.names = T, sep = '\t', append = TRUE)
close(filecon)
}
