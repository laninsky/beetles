library(data.table)

linenos <- as.matrix(fread("linestosample",integer64="numeric"))

print("I have read in the line numbers of the sequence titles")

newmatrix <- matrix("NA",ncol=1,nrow=(dim(linenos)[1]*4))
newmatrix[(seq(1,dim(newmatrix)[1],4)),1] <- linenos[,1]
newmatrix[(seq(2,dim(newmatrix)[1],4)),1] <- linenos[,1]+1
newmatrix[(seq(3,dim(newmatrix)[1],4)),1] <- linenos[,1]+2
newmatrix[(seq(4,dim(newmatrix)[1],4)),1] <- linenos[,1]+3

rm(linenos)

print("I have calculated the line numbers that I need to subsample from your fastq file")

filename <- as.matrix(read.table("filename"))

fastqfile <- as.matrix(fread(filename,integer64="character",sep="\n"))

print("I have read in your fastq file")

output <- fastqfile[as.numeric(newmatrix[,1]),]

write.table(output, "downsampled.fastq",quote=FALSE, col.names=FALSE,row.names=FALSE,append=TRUE)

print("I have downsampled your fastq file and written it as downsampled.fastq")
