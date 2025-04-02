args = commandArgs(trailingOnly=TRUE)


barcodes<-read.table(args[1],h=T,sep="\t")

barcodeStacks<-data.frame(sample=barcodes$Sample,
                          barcode1=paste0(barcodes$Barcode_R1   ,"C"),
                          barcode2=paste0(barcodes$Barcode_R2,"C"))
write.table(barcodeStacks,paste0(args[2]),row.names = F,col.names = F,quote=F,sep="\t")
