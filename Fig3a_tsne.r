rm(list=ls())

library(tsne)
library('tidyverse')
library(Rtsne)


work.dir = "~/source code/Data_Fig3a"
setwd(work.dir)

cluster_file = "Input_data_fig3a.txt"     


cluster.data = read.delim(cluster_file,stringsAsFactors = F,quote="",sep="\t")
cluster.matrix = as.matrix(cluster.data[4:(ncol(cluster.data))])


cluster.scaled.mat = scale(cluster.matrix, center=T, scale = TRUE)
Type = cluster.data$type
Cancer.type = cluster.data$cancer
Sample.type = cluster.data$sample_type
cluster.tsne = tsne(cluster.scaled.mat,k=2,perplexity = 30)  #pan_cancer 


cluster.tsne.F = as.data.frame(cluster.tsne)
cluster.tsne.F = cbind(cluster.tsne.F,Type,Cancer.type,Sample.type)

write.table(cluster.tsne.F,paste(cluster_file,"_xy.txt"),row.names = F, quote = F,sep="\t")

