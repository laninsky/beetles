removing_phylip_taxon <- function(working_dir,taxon_name) {
  cat("This function takes a folder of phylip files and writes out new versions of them removing a defined taxon")
  cat("Run this function by: removing_phylip_taxon(working_dir,taxon_name)\n")
  cat("Where working_dir = the path to a folder containing your phylip files and taxon_name the name of the taxon you wish to remove from the phylip files\n")
  cat('e.g. removing_phylip_taxon("/Users/alanaalexander/Dropbox/UCE_data_for_Alana/Spanglerogyrus_removal/50perc_phylip_alignments_w_Spangy","sle1702_sp")\n')
  cat("Your new phylip files will be written out to a folder in your directory called taxon_name_removed\n")
  cat("e.g. sle1702_sp_removed\n")
  
  # getting a list of the files from the folder which match a phylip file
  file_list <- list.files(working_dir, pattern=".phy*")
  
  # creating output directory
  dir.create(paste(taxon_name,"_removed",sep=""))

  # For each of the phylip files...
  for (i in file_list) {
    
    # Reading in the file
    tempfile <- readLines(file.path(working_dir,i))
    
    # Getting the total number of taxa in the file
    no_taxa <- unlist(strsplit(tempfile[1], " "))
    no_taxa <- as.numeric(no_taxa[length(no_taxa)-1])
    
    # Getting rid of any blank linkes
    tempfile <- tempfile[which(grepl("^$",tempfile)==FALSE)]
    
    # Converting the file to non-interleaved/unwrapped if interleaved
    if (length(tempfile)>(no_taxa+1)) {
      # Calculating the number of times each sequence has been wrapped
      wrapped <- (length(tempfile)-1)/no_taxa
      
      # Then for each taxon, pasting its wrapped lines to the first instance it occurs 
      for (j in 2:(no_taxa+1)) {
        for (k in 1:(wrapped-1)) {
          tempfile[j] <- paste(tempfile[j],tempfile[k*no_taxa+j],sep="")
        }
        # Removing spaces in the sequence, and assuming phylip 10-character name limit
        tempjname <- unlist(strsplit(tempfile[j]," "))[1]
        tempjseq <- unlist(strsplit(tempfile[j]," "))[-1]
        tempjseq <- paste(tempjseq[which(tempjseq!="")],collapse="")
        tempfile[j] <- paste(tempjname,paste(rep(" ",(11-nchar(tempjname))),collapse=""),tempjseq,sep="")
      }
      # Taking just the first sequencing line per taxon as this now has all sequencing data for it
      tempfile <- tempfile[1:(no_taxa+1)]
    }  
    
    # Removing taxon
    tempfile[-which(grepl(taxon_name,tempfile))]
    
    # Adjusting number of taxa in first line
    tempfile[1] <- gsub(paste(no_taxa," ",sep=""),paste(no_taxa-1," ",sep=""),tempfile[1])

    # Writing out file
    write.table(tempfile,file.path(paste(taxon_name,"_removed",sep=""),i),quote=FALSE,row.names = FALSE,col.names = FALSE)
    
  }
}  
