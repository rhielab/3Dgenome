library(ggplot2)
library(dplyr)
library(ggpubr)


enhancer_enhancer <- read.table("Enhancer_in_loop_with_enhancer.bed")
enhancer_enhancer['Loop_Category'] <- "Enhancer-Enhancer"
enhancer_promoter <- read.table("Enhancer_in_loop_with_promoter.bed")
enhancer_promoter['Loop_Category'] <- "Enhancer-Promoter"
enhancer_ndr <- read.table("Enhancer_in_loop_with_ndr.bed")
enhancer_ndr['Loop_Category'] <- "Enhancer-NDR"
enhancer_insulator <- read.table("Enhancer_in_loop_with_insulator.bed")
enhancer_insulator['Loop_Category'] <- "Enhancer-Insulator"
enhancer_h3k27me3 <- read.table("Enhancer_in_loop_with_h3k27me3.bed")
enhancer_h3k27me3['Loop_Category'] <- "Enhancer-H3K27me3"
enhancer_h3k9me3 <- read.table("Enhancer_in_loop_with_h3k9me3.bed")
enhancer_h3k9me3['Loop_Category'] <- "Enhancer-H3K9me3"
enhancer_none <- read.table("Enhancer_in_loop_with_none.bed")
enhancer_none['Loop_Category'] <- "Enhancer-None"



all_enhancer <- rbind(enhancer_enhancer,enhancer_promoter, enhancer_ndr, enhancer_insulator, enhancer_h3k27me3, enhancer_h3k9me3, enhancer_none )




jpeg("Enhancer_itself_by_loop_category.jpg",, width = 1600, height = 1600)
ggplot(data = all_enhancer, aes(y=V7, x=Loop_Category)) + geom_violin( fill = "grey") + ylim(0,100) + stat_summary(fun = "mean",geom = "point",color = "red", size = 5) + stat_summary(fun = "median", geom = "point",color = "blue", size = 5)+ ylab("Signal Value") + xlab("Enhancer Categories") + theme(axis.text =  element_text(size = 5)) + scale_fill_grey() + theme_bw() + theme(  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
dev.off()

insulator_insulator <- read.table("Insulator_in_loop_with_insulator.bed")
insulator_insulator['Loop_Category'] <- "Insulator-Insulator"
insulator_promoter <- read.table("Insulator_in_loop_with_promoter.bed")
insulator_promoter['Loop_Category'] <- "Insulator-Promoter"
insulator_ndr <- read.table("Insulator_in_loop_with_ndr.bed")
insulator_ndr['Loop_Category'] <- "Insulator-NDR"
insulator_enhancer <- read.table("Insulator_in_loop_with_enhancer.bed")
insulator_enhancer['Loop_Category'] <- "Insulator-Enhancer"
insulator_h3k27me3 <- read.table("Insulator_in_loop_with_h3k27me3.bed")
insulator_h3k27me3['Loop_Category'] <- "Insulator-H3K27me3"
insulator_h3k9me3 <- read.table("Insulator_in_loop_with_h3k9me3.bed")
insulator_h3k9me3['Loop_Category'] <- "Insulator-H3K9me3"
insulator_none <- read.table("Insulator_in_loop_with_none.bed")
insulator_none['Loop_Category'] <- "Insulator-None"



all_insulator <- rbind(insulator_insulator,insulator_promoter, insulator_ndr, insulator_enhancer, insulator_h3k27me3, insulator_h3k9me3, insulator_none )




jpeg("Insulator_itself_by_loop_category.jpg",, width = 1600, height = 1600)
ggplot(data = all_insulator, aes(y=V7, x=Loop_Category)) + geom_violin( fill = "grey") + ylim(0,100)+ stat_summary(fun = "mean",geom = "point",color = "red", size = 5) + stat_summary(fun = "median", geom = "point",color = "blue", size = 5) + ylab("Signal Value") + xlab("Insulator Categories") + theme(axis.text =  element_text(size = 5)) + scale_fill_grey() + theme_bw() + theme(  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
dev.off()


ndr_promoter <- read.table("NDR_in_loop_with_promoter.bed")
ndr_promoter['Loop_Category'] <- "NDR-Promoter"
ndr_insulator <- read.table("NDR_in_loop_with_insulator.bed")
ndr_insulator['Loop_Category'] <- "NDR-Insulator"
ndr_enhancer <- read.table("NDR_in_loop_with_enhancer.bed")
ndr_enhancer['Loop_Category'] <- "NDR-Enhancer"
ndr_h3k27me3 <- read.table("NDR_in_loop_with_h3k27me3.bed")
ndr_h3k27me3['Loop_Category'] <- "NDR-H3K27me3"
ndr_h3k9me3 <- read.table("NDR_in_loop_with_h3k9me3.bed")
ndr_h3k9me3['Loop_Category'] <- "NDR-H3K9me3"
ndr_none <- read.table("NDR_in_loop_with_none.bed")
ndr_none['Loop_Category'] <- "NDR-None"
ndr_ndr <- read.table("NDR_in_loop_with_ndr.bed")
ndr_ndr['Loop_Category'] <- "NDR-NDR"



all_ndr <- rbind(ndr_ndr,ndr_promoter, ndr_insulator, ndr_enhancer, ndr_h3k27me3, ndr_h3k9me3, ndr_none )




jpeg("NDR_itself_by_loop_category.jpg", width = 1600, height = 1600)
ggplot(data = all_ndr, aes(y=V5, x=Loop_Category)) + geom_violin( fill = "grey") + ylim(2.5,10) + ylab("Signal Value") + stat_summary(fun = "mean",geom = "point",color = "red", size = 5) + stat_summary(fun = "median", geom = "point",color = "blue", size = 5)+ xlab("NDR Categories") + theme(axis.text =  element_text(size = 5)) + scale_fill_grey() + theme_bw() + theme(  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
dev.off()

h3k27me3_promoter <- read.table("H3K27me3_in_loop_with_promoter.bed")
h3k27me3_promoter['Loop_Category'] <- "H3K27me3-Promoter"
h3k27me3_insulator <- read.table("H3K27me3_in_loop_with_insulator.bed")
h3k27me3_insulator['Loop_Category'] <- "H3K27me3-Insulator"
h3k27me3_enhancer <- read.table("H3K27me3_in_loop_with_enhancer.bed")
h3k27me3_enhancer['Loop_Category'] <- "H3K27me3-Enhancer"
h3k27me3_ndr <- read.table("H3K27me3_in_loop_with_NDR.bed")
h3k27me3_ndr['Loop_Category'] <- "H3K27me3-NDR"
h3k27me3_h3k9me3 <- read.table("H3K27me3_in_loop_with_h3k9me3.bed")
h3k27me3_h3k9me3['Loop_Category'] <- "H3K27me3-H3K9me3"
h3k27me3_none <- read.table("H3K27me3_in_loop_with_none.bed")
h3k27me3_none['Loop_Category'] <- "H3K27me3-None"
h3k27me3_h3k27me3 <- read.table("H3K27me3_in_loop_with_h3k27me3.bed")
h3k27me3_h3k27me3['Loop_Category'] <- "H3K27me3-H3K27me3"


all_h3k27me3 <- rbind(h3k27me3_h3k27me3, h3k27me3_promoter, h3k27me3_insulator, h3k27me3_enhancer, h3k27me3_ndr, h3k27me3_h3k9me3, h3k27me3_none )




jpeg("H3K27me3_itself_by_loop_category.jpg", width = 1600, height = 1600)
ggplot(data = all_h3k27me3, aes(y=V7, x=Loop_Category)) + geom_violin( fill = "grey") + ylim(1.5,3) + ylab("Signal Value") + stat_summary(fun = "mean",geom = "point",color = "red", size = 5) + stat_summary(fun = "median", geom = "point",color = "blue", size = 5)+ xlab("H3K27me3 Categories") + theme(axis.text =  element_text(size = 5)) + scale_fill_grey() + theme_bw() + theme(  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
dev.off()
h3k9me3_promoter <- read.table("H3K9me3_in_loop_with_promoter.bed")
h3k9me3_promoter['Loop_Category'] <- "H3K9me3-Promoter"
h3k9me3_insulator <- read.table("H3K9me3_in_loop_with_insulator.bed")
h3k9me3_insulator['Loop_Category'] <- "H3K9me3-Insulator"
h3k9me3_enhancer <- read.table("H3K9me3_in_loop_with_enhancer.bed")
h3k9me3_enhancer['Loop_Category'] <- "H3K9me3-Enhancer"
h3k9me3_ndr <- read.table("H3K9me3_in_loop_with_NDR.bed")
h3k9me3_ndr['Loop_Category'] <- "H3K9me3-NDR"
h3k9me3_h3k27me3 <- read.table("H3K9me3_in_loop_with_h3k27me3.bed")
h3k9me3_h3k27me3['Loop_Category'] <- "H3K9me3-H3K27me3"
h3k9me3_none <- read.table("H3K9me3_in_loop_with_none.bed")
h3k9me3_none['Loop_Category'] <- "H3K9me3-None"
h3k9me3_h3k9me3 <- read.table("H3K9me3_in_loop_with_h3k9me3.bed")
h3k9me3_h3k9me3['Loop_Category'] <- "H3K9me3-H3K9me3"


all_h3k9me3 <- rbind(h3k9me3_h3k9me3,h3k9me3_promoter, h3k9me3_insulator, h3k9me3_enhancer, h3k9me3_ndr, h3k9me3_h3k27me3, h3k9me3_none )




jpeg("H3K9me3_itself_by_loop_category.jpg", width = 1600, height = 1600)
ggplot(data = all_h3k9me3, aes(y=V7, x=Loop_Category)) + geom_violin( fill = "grey") + ylim(1,4) + ylab("Signal Value") + stat_summary(fun = "mean",geom = "point",color = "red", size = 5) + stat_summary(fun = "median", geom = "point",color = "blue", size = 5)+ xlab("H3K9me3 Categories") + theme(axis.text =  element_text(size = 5)) + scale_fill_grey() + theme_bw() + theme(  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
dev.off()

