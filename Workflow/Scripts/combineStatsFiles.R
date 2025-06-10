rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
#directory<-"/vol/ribesecology/nielsw/constant/msGBS/Output/Analysis/Bowtie/perSample"
#output<-"/vol/ribesecology/nielsw/constant/msGBS/Output/Analysis/Bowtie/stats.tsv"
directory<-args[1]
output<-args[2]


library("purrr")
library(tidyverse)

fileList<-list.files(directory,pattern="tsv",full.names = T)

csvList<-lapply(fileList,read.table,h=T)
csvListSub<-csvList[map(csvList, function(x) dim(x)[1]) > 0]
removedSamples<-fileList[map(csvList, function(x) dim(x)[1]) == 0]
removedSamples<-sub(".tsv","",gsub(".*/","",removedSamples))
out<-csvListSub %>%
  reduce(full_join, by = c("Species","Locus"))
out<-out[order(out$Locus),]
out[is.na(out)]<-0
colnames(out)<-sub("^X","",colnames(out))

print(removedSamples)
write.table(out,output,row.names=F, sep='\t', quote=F)
write.table(removedSamples,paste(output,"removedSamples.txt"))

