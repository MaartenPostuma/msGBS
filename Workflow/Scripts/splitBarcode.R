args = commandArgs(trailingOnly=TRUE)


barcodes<-read.table(args[1],h=T,sep="\t")

barcodeStacks<-data.frame(barcode1=paste0(barcodes$Barcode_R1   ,"C"),
                          barcode2=paste0(barcodes$Barcode_R2,"C"),
                          sample=barcodes$Sample)
write.table(barcodeStacks,paste0(args[2]),row.names = F,col.names = F,quote=F,sep="\t")
