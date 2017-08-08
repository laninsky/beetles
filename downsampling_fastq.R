library(data.table)

linenos <- as.matrix(fread("linestosample",integer64="numeric"))

newmatrix <- matrix("NA",ncol=1,nrow=(dim(linenos)[1]*4))
newmatrix[(seq(1,dim(newmatrix)[1],4)),1] <- linenos[,1]
newmatrix[(seq(2,dim(newmatrix)[1],4)),1] <- linenos[,1]+1
newmatrix[(seq(3,dim(newmatrix)[1],4)),1] <- linenos[,1]+2
newmatrix[(seq(4,dim(newmatrix)[1],4)),1] <- linenos[,1]+3

rm(linenos)

filename <- as.matrix(read.table("filename"))

fastqfile <- as.matrix(fread(filename,integer64="character"))




output <- fastqfile[as.numeric(newmatrix[,1]),]

write.table(output, "downsampled.fastq",quote=FALSE, col.names=FALSE,row.names=FALSE)
