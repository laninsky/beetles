library(data.table)

linenos <- as.matrix(fread("linestosample",integer64="numeric"))

newmatrix <- matrix("NA",ncol=1,nrow=(dim(linenos)[1]*4))
newmatrix[(seq(1,dim(newmatrix)[1],4)),1] <- linenos[,1]
newmatrix[(seq(2,dim(newmatrix)[1],4)),1] <- linenos[,1]+1

#check this coding is correct
newmatrix[(seq(3,dim(newmatrix)[1],4)),1] <- linenos[,1]+2
#check this coding is correct
newmatrix[(seq(4,dim(newmatrix)[1],4)),1] <- linenos[,1]+3

rm(linenos)

file <- as.matrix(fread("filename"))

output <- file[as.numeric(newmatrix[,1]),]

write.table(output, "downsampled.fastq",quote=FALSE, col.names=FALSE,row.names=FALSE)
