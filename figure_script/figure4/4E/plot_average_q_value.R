library(dplyr)
library(fuzzyjoin)
library(ggplot2)

chr_limit <-  c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")

ee_loop <- read.table("OutputPrefix_Enhancer_Enhancer_Loop.tsv")
eh27_loop <- read.table("OutputPrefix_Enhancer_H3K27me3_Loop.tsv")
eh9_loop <- read.table("OutputPrefix_Enhancer_H3K9me3_Loop.tsv")
ei_loop <- read.table("OutputPrefix_Enhancer_Insulator_Loop.tsv")
en_loop <- read.table("OutputPrefix_Enhancer_NDR_Loop.tsv")
ex_loop <- read.table("OutputPrefix_Enhancer_None_Loop.tsv")
h27h27_loop <- read.table("OutputPrefix_H3K27me3_H3K27me3_Loop.tsv")
h27h9_loop <- read.table("OutputPrefix_H3K27me3_H3K9me3_Loop.tsv")
h27x_loop <- read.table("OutputPrefix_H3K27me3_None_Loop.tsv")
h9h9_loop <- read.table("OutputPrefix_H3K9me3_H3K9me3_Loop.tsv")
h9x_loop <- read.table("OutputPrefix_H3K9me3_None_Loop.tsv")
ih27_loop <- read.table("OutputPrefix_Insulator_H3K27me3_Loop.tsv")
ih9_loop <- read.table("OutputPrefix_Insulator_H3K9me3_Loop.tsv")
ii_loop <- read.table("OutputPrefix_Insulator_Insulator_Loop.tsv")
in_loop <- read.table("OutputPrefix_Insulator_NDR_Loop.tsv")
ix_loop <- read.table("OutputPrefix_Insulator_None_Loop.tsv")
nh27_loop <- read.table("OutputPrefix_NDR_H3K27me3_Loop.tsv")
nh9_loop <- read.table("OutputPrefix_NDR_H3K9me3_Loop.tsv")
nn_loop <- read.table("OutputPrefix_NDR_NDR_Loop.tsv")
nx_loop <- read.table("OutputPrefix_NDR_None_Loop.tsv")
xx_loop <- read.table("OutputPrefix_None_None_Loop.tsv")
pe_loop <- read.table("OutputPrefix_Promoter_Enhancer_Loop.tsv")
ph27_loop <- read.table("OutputPrefix_Promoter_H3K27me3_Loop.tsv")
ph9_loop <- read.table("OutputPrefix_Promoter_H3K9me3_Loop.tsv")
pi_loop <- read.table("OutputPrefix_Promoter_Insulator_Loop.tsv")
pn_loop <- read.table("OutputPrefix_Promoter_NDR_Loop.tsv")
px_loop <- read.table("OutputPrefix_Promoter_None_Loop.tsv")
pp_loop <- read.table("OutputPrefix_Promoter_Promoter_Loop.tsv")

ee_loop['Loop_category'] <- "Enhancer\nEnhancer"
eh27_loop['Loop_category'] <- "Enhancer\nH3K27me3"
eh9_loop['Loop_category'] <- "Enhancer\nH3K9me3"
ei_loop['Loop_category'] <- "Enhancer\nInsulator"
en_loop['Loop_category'] <- "Enhancer\nNDR"
ex_loop['Loop_category'] <- "Enhancer\nNone"
h27h27_loop['Loop_category'] <- "H3K27me3\nH3K27me3"
h27h9_loop['Loop_category'] <- "H3K27me3\nH3K9me3"
h27x_loop['Loop_category'] <- "H3K27me3\nNone"
h9h9_loop['Loop_category'] <- "H3K9me3\nH3K9me3"
h9x_loop['Loop_category'] <- "H3K9me3\nNone"
ih27_loop['Loop_category'] <- "Insulator\nH3K27me3"
ih9_loop['Loop_category'] <- "Insulator\nH3K9me3"
ii_loop['Loop_category'] <- "Insulator\nInsulator"
in_loop['Loop_category'] <- "Insulator\nNDR"
ix_loop['Loop_category'] <- "Insulator\nNone"
nh27_loop['Loop_category'] <- "NDR\nH3K27me3"
nh9_loop['Loop_category'] <- "NDR\nH3K9me3"
nn_loop['Loop_category'] <- "NDR\nNDR"
nx_loop['Loop_category'] <- "NDR\nNone"
xx_loop['Loop_category'] <- "None\nNone"
pe_loop['Loop_category'] <- "Promoter\nEnhancer"
ph27_loop['Loop_category'] <- "Promoter\nH3K27me3"
ph9_loop['Loop_category'] <- "Promoter\nH3K9me3"
pi_loop['Loop_category'] <- "Promoter\nInsulator"
pn_loop['Loop_category'] <- "Promoter\nNDR"
px_loop['Loop_category'] <- "Promoter\nNone"
pp_loop['Loop_category'] <- "Promoter\nPromoter"


all_of_interaction <- rbind( ei_loop , ih27_loop, ii_loop, pe_loop, pi_loop )

jpeg("OutputPrefix_Loop_Category_Top5_Q_Value_Violin_Plot.jpg", width = 1600, height = 1600)
ggplot(data = all_of_interaction, aes(y=V7, x =reorder(Loop_category,V7,FUN = median))) + geom_violin( fill = "grey") + stat_summary(fun = "mean",geom = "point",color = "red", size = 5) + stat_summary(fun = "median", geom = "point",color = "blue", size = 5)+ ylab("MicroC_3_Billion_5kb Interaction Distribution") + xlab("Loop Category") + ylim(0,0.20) + theme(axis.text =  element_text(size = 10)) + scale_fill_grey() + theme_bw() + theme(  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())

hic_anova <- aov(V7 ~ Loop_category, data = all_of_interaction)
summary(hic_anova)

t.test(ii_loop$V7, ih27_loop$V7)
