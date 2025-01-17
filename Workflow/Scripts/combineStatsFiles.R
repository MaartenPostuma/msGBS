rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
directory<-args[1]
output<-args[2]
library("purrr")
library(tidyverse)

fileList<-list.files(directory,pattern="csv",full.names = T)

csvList<-lapply(fileList,read.table,h=T)
csvListSub<-csvList[map(csvList, function(x) dim(x)[1]) > 0]
removedSamples<-fileList[map(csvList, function(x) dim(x)[1]) == 0]
removedSamples<-sub(".csv","",gsub(".*/","",removedSamples))
out<-csvListSub %>%
  reduce(full_join, by = c("Species","Locus"))
out<-out[order(out$Locus),]
out[is.na(out)]<-0
colnames(out)<-sub("^X","",colnames(out))

print(removedSamples)
write.csv(out,output,row.names=F)
write.csv(removedSamples,paste(output,"removedSamples.txt"))

