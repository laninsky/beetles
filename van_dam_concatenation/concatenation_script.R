# Assumptions: taxa are named the same way in all of the uce files, 
# and this matches the way they are named in uce_summary_file_path
# all taxa are present in all phylip files pointed to in phylip_w_missing_path
# (e.g. they are phylip files w missing data added). If there are
# files in uce_summary_file_path that aren't present in uce_summary_file_path
# the code should handle this OK

# Pathway to the *uce_type_summary.txt file
#uce_summary_file_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/All_beetle_probeset_of_Faircloth/Coleoptera_Tricas_uce_type_summary.txt"
uce_summary_file_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/Adephaga_dataset/Adephaga_Tricas_uce_type_summary.txt"

# Pathway to the phylip input files
#phylip_w_missing_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/All_beetle_probeset_of_Faircloth/50perc_internal_w_missing_phylip"
phylip_w_missing_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/Adephaga_dataset/50perc_internal_w_missing_phylip"

# Where you would like the outfiles written to
#outfile_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/All_beetle_probeset_of_Faircloth/van_dam_concatenated"
outfile_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/Adephaga_dataset/van_dam_concatenated"

# Pathway to the phylip files without missing data
#phylip_without_missing_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/All_beetle_probeset_of_Faircloth/50perc_internal_phylip"
phylip_without_missing_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/Adephaga_dataset/50perc_internal_phylip"

library(tidyverse)

# Creating the output directory
dir.create(file.path(outfile_path), showWarnings = FALSE)
dir.create(file.path(outfile_path,"non_concat_loci"), showWarnings = FALSE)

# Reading in the uce summary file
uce_summary <- read_delim(uce_summary_file_path,delim="\t",col_names = FALSE)

# Extracting the gene name from the extended description
# Simplifying to just the uce name and gene_name
uce_summary <- uce_summary %>% 
  mutate(gene_name=gsub("gene(ID=gene-","",gsub(").*","",X6),fixed=TRUE)) %>% 
  select(X1,gene_name)

loci_to_not_copy <- NULL

for (i in unique(uce_summary$gene_name)) {
  # Extracting the UCE names for the gene name we are up to
  uce_subsets <- uce_summary$X1[which(uce_summary$gene_name==i)]
  
  # Getting a variable together to concatenate the loci across
  out_uce <- NULL
  
  # For each of the remaining UCEs that map to that gene
  for (j in uce_subsets) {
    
    # Matching the uce name that is in the folder
    uce_name <- list.files(phylip_w_missing_path, pattern=paste(j,"\\.",sep=""))[1]
    
    if (!is.na(uce_name)) {
      
      # Loci that are written out here as part of a gene do not need to be copied across
      loci_to_not_copy <- c(loci_to_not_copy,uce_name)
      
      # Reading in this UCE and removing all empty lines
      uce <- readLines(paste(phylip_w_missing_path,"/",uce_name,sep=""))
      uce <- uce[-which(uce=="")]
      
      # Number of taxa
      numtaxa <- as.numeric(gsub(" .*","",gsub("^ *","",uce[1])))
      
      # Extracting the taxa names
      uce_taxa <-  gsub(" .*","",uce[2:(numtaxa+1)])
      
      # Getting the first sequence line together
      seq_line <-  uce[2:(numtaxa+1)]
      for (k in 1:length(seq_line)) {
        seq_line[k] <- gsub(" ","",gsub(uce_taxa[k],"",seq_line[k]))
      }
      
      # Concatenating the rest of the sequences
      for (k in 1:length(seq_line)) {
        seq_line[k] <- gsub(" ","",paste(seq_line[k],paste(gsub(" ","",uce[seq((numtaxa+1+k),length(uce),numtaxa)]),collapse=""),collapse=""))
      }
      
      # If out_uce doesn't have any data in it yet
      if (is.null(out_uce)) {
        # Writing out this locus to out_uce
        out_uce <- as_tibble(cbind(uce_taxa,seq_line)) %>% arrange(uce_taxa)
      # If out_uce already has something in it, then we need to concatenate this current loci to it  
      } else {
        # Ordering this current loci the same way as in out_uce
        temp_out_uce <- as_tibble(cbind(uce_taxa,seq_line)) %>% arrange(uce_taxa)
        
        # Going through and adding sequence from each loci
        for (k in 1:dim(out_uce)[1]) {
          out_uce$seq_line[k] <- gsub(" ","",paste(out_uce$seq_line[k],temp_out_uce$seq_line[k],collapse = ""))
        }
      }
    }
  }
  # What to do to write out this specific locus
  if (!(is.null(out_uce))) {
    
    # Filtering out taxa that are entirely composed of missing data
    out_uce <- out_uce %>% mutate(missing_data=nchar(gsub("?","",seq_line,fixed=TRUE))) %>% filter(missing_data!=0)
    
    # Recalculating the number of missing taxa
    numtaxa <- dim(out_uce)[1]
    
    # Getting the sequence length
    seq_length <- nchar(out_uce$seq_line[1])
    
    # Calculating the number of characters to separate uce_taxa and seq_line
    num_spaces <- max(nchar(out_uce$uce_taxa))+1
    
    # Pasting the names and sequence together
    to_write <- out_uce %>% 
      rowwise() %>% 
      mutate(spaces=paste(rep(" ",num_spaces-nchar(uce_taxa)),collapse="")) %>% 
      mutate(final_out=paste(uce_taxa,spaces,seq_line,collapse="")) %>% select(final_out)
    
    # Adding the information line to the top
    to_write <- rbind(paste(" ",numtaxa," ",seq_length,collape=""),as.matrix(to_write$final_out))
    
    # Writing out the data
    write.table(file=gsub(" ","",paste(outfile_path,"/",i,".phylip",collapse="")),x=to_write,quote=FALSE,col.names = FALSE,row.names = FALSE)
  }
}

# Writing out the list of files that has been concatenated
write.table(file=gsub(" ","",paste(outfile_path,"/","loci_that_were_concatenated.txt",collapse="")),x=loci_to_not_copy,quote=FALSE,col.names = FALSE,row.names = FALSE)

# Getting the list of loci to copy into the concatenated folder that AREN'T already represented in
# the concatenated loci
copy_list <- list.files(phylip_without_missing_path)[!(list.files(phylip_without_missing_path) %in% loci_to_not_copy)]
# Excluding anything that isn't a phylip
copy_list <- copy_list[grep(pattern = "phylip",copy_list)]

for (i in copy_list) {
  file.copy(from = gsub(" ","",paste(phylip_without_missing_path,"/",i,collapse="")),to=gsub(" ","",paste(outfile_path,"/non_concat_loci",collapse="")))
}

