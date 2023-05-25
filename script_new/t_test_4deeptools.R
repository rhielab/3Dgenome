#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Two output files from deeptools should be used!", call.=FALSE)
} else {
  data1 <- read.table(file = args[1], skip = 1, sep = '\t')
  data2 <- read.table(file = args[2], skip = 1, sep = '\t')
}

print(t.test(data1[96:116], data2[96:116], na.rm = TRUE))
