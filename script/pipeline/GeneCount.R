library(dplyr)
expdesign <- read.csv(file = '/scratch/beoungle/probing_rnaseq/c42b/output/c42b/c42b_MIX_expdesign.csv', header = T, check.names = F, row.names = 1)
quantification <- 'HTseq'

if(quantification == 'HTseq') {
a <- grep('HTseq',list.dirs('/scratch/beoungle/probing_rnaseq/c42b/output/c42b', recursive = F), value = T) %>% list.files(pattern = '.counts')
DIR <- grep('HTseq',list.dirs('/scratch/beoungle/probing_rnaseq/c42b/output/c42b', recursive = F), value = T) %>% list.files(pattern = '.counts', full.names = T)
final <- read.table(file = DIR[1], header = F, sep = '\t', row.names = 1)
for (i in 2:length(DIR)) {
	new.data <- read.table(file = DIR[i], header = F, sep = '\t', row.names = 1)
    final <- cbind(final, new.data)
}
names(final) <- a
}

if(quantification == 'STAR') {
a <- grep('STAR',list.dirs('/scratch/beoungle/probing_rnaseq/c42b/output/c42b', recursive = F), value = T) %>% list.files(pattern = 'ReadsPerGene.out.tab')
DIR <- grep('STAR',list.dirs('/scratch/beoungle/probing_rnaseq/c42b/output/c42b', recursive = F), value = T) %>% list.files(pattern = 'ReadsPerGene.out.tab', full.names = T)
final <- read.table(file = DIR[1], header = F, sep = '\t', row.names = 1)
final <- final[-c(1:4),]
final <- subset(final, select = 1)
for (i in 2:length(DIR)) {
    new.data <- read.table(file = DIR[i], header = F, check.names = F, sep = '\t', row.names = 1)
    new.data <- new.data[-c(1:4),]
    new.data <- subset(new.data, select = 1)
    final <- cbind(final, new.data)
}
names(final) <- a
}

if(quantification == 'Rsem') {
a <- grep('Rsem',list.dirs('/scratch/beoungle/probing_rnaseq/c42b/output/c42b', recursive = F), value = T) %>% list.files(pattern = 'genes.results')
DIR <- grep('Rsem',list.dirs('/scratch/beoungle/probing_rnaseq/c42b/output/c42b', recursive = F), value = T) %>% list.files(pattern = 'genes.resultsb', full.names = T)
final <- read.table(file = DIR[1], header = T, check.names = F, sep = '\t', row.names = 1)
final <- subset(final, select = 4)
n <- length(a)
for (i in 2:n) {
new.data <- read.table(file = DIR[i], header = T, check.names = F, sep = '\t', row.names = 1)
new.data <- subset(new.data, select = 4)
final <- cbind(final,new.data)
}
names(final) <- a
coluname <- colnames(final)
rowname <- rownames(final)
final <- apply(final, 2, as.integer)
rownames(final) <- rowname
colnames(final) <- coluname
final <- as.data.frame(final)
}

filename <- rownames(expdesign) 
element <- strsplit(filename[1],'[.]') %>% unlist()
if(element[length(element)] == 'gz' && element[length(element)-1] == 'fq') {
extension = '.fq.gz';nn = 2
} else if(element[length(element)] == 'gz' && element[length(element)-1] == 'fastq') {
extension = '.fastq.gz';nn = 2
} else if(element[length(element)] == 'fq') {
extension = '.fq'; nn = 1
} else if(element[length(element)] == 'fastq') {
extension = '.fastq'; nn = 1
}
Prefix <- lapply(filename, function(x){element2 <- strsplit(x,'[.]') %>% unlist(); 
paste(element2[1:(length(element2)-nn) ], collapse = '.')}) %>% unlist()
Prefix <- sort(Prefix)
filename <- sort(filename)
names(final) <- filename[lapply(Prefix, function(x){grep(x,names(final))} ) %>% unlist()]

write.csv(final, file = '/scratch/beoungle/probing_rnaseq/c42b/output/c42b/HTseq-GenesCountsTable.csv', quote = F, row.names = T)