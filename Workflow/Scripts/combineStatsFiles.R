rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
#args<-c()
#args[1]<-"/vol/ribesecology/nielsw/constant/constant_2023-2024/output/Analysis/Bowtie/monos/perSample/"
#args[2]<-"/vol/ribesecology/nielsw/constant/constant_2023-2024/output/Analysis/Bowtie/stats.tsv"
#args[3]<-"/vol/ribesecology/nielsw/constant/constant_2023-2024/output/Analysis/Bowtie/nonmonos/perSample/"

if (grepl("perSample",args[3])){
directoryMono<-args[1]
output<-args[2]
directoryNonMono<-args[3]
fileListMono<-list.files(directoryMono,pattern="tsv",full.names = T)
directoryNonMono<-list.files(directoryNonMono,pattern="tsv",full.names = T)
fileList<-c(fileListMono,directoryNonMono)
 }else{
directoryMono<-args[1]
output<-args[2]
fileList<-list.files(directoryMono,pattern="tsv",full.names = T)
}

library("purrr")
library(tidyverse)
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

