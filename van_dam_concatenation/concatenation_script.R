uce_summary_file_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/Adephaga_dataset/Adephaga_Tricas_uce_type_summary.txt"

phylip_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/Adephaga_dataset/50perc_internal_phylip"

library(tidyverse)

# Reading in the uce summary file
uce_summary <- read_delim(uce_summary_file_path,delim="\t",col_names = FALSE)

# Extracting the gene name from the extended description
# Simplifying to just the uce name and gene_name
uce_summary <- uce_summary %>% 
  mutate(gene_name=gsub("gene(ID=gene-","",gsub(").*","",X6),fixed=TRUE)) %>% 
  select(X1,gene_name)

for (i in unique(uce_summary$gene_name)) {
  # Extracting the UCE names for the gene name we are up to
  uce_subsets <- uce_summary$X1[which(uce_summary$gene_name==i)]
  
  # Getting the name of the first uce with this specific gene name
  first_uce_name <- list.files(phylip_path, pattern=paste(uce_subsets[1],"\\.",sep=""))[1]
  
  # Reading in this UCE and removing all empty lines
  first_uce <- readLines(paste(phylip_path,"/",first_uce_name,sep=""))
  first_uce <- first_uce[-which(first_uce=="")]
  
  # Number of taxa
  numtaxa <- as.numeric(gsub(" .*","",gsub("^ *","",first_uce[1])))
  
  # Extracting the taxa names
  first_uce_taxa <-  gsub(" .*","",first_uce[2:(numtaxa+1)])
  
  # Getting the first sequence line together
  seq_line <-  first_uce[2:(numtaxa+1)]
  for (j in 1:length(seq_line)) {
    seq_line[j] <- gsub(" ","",gsub(first_uce_taxa[j],"",seq_line[j]))
  }
  
  # Concatenating the rest of the sequences
  for (j in 1:length(seq_line)) {
    seq_line[j] <- gsub(" ","",paste(seq_line[j],paste(gsub(" ","",first_uce[seq((numtaxa+1+j),length(first_uce),numtaxa)]),collapse=""),collapse=""))
  }
  
  out_uce <- as_tibble(cbind(first_uce_taxa,seq_line)) %>% arrange(first_uce_taxa)
                         
  
   UP TO HEREEEEEEEEEEE
  
  
  # For each of the remaining UCEs that map to that gene
  for (j in uce_subsets[-1]) {
    # Getting the name of the this uce
    uce_name <- list.files(w_missing_phylip_path, pattern=paste(j,"\\.",sep=""))
    
    # Reading in this UCE and removing all empty lines
    first_uce <- readLines(paste(w_missing_phylip_path,"/",first_uce_name,sep=""))
    first_uce <- first_uce[-which(first_uce=="")]
    
    # Number of taxa
    numtaxa <- as.numeric(gsub(" .*","",gsub("^ *","",first_uce[1])))
    
    # Extracting the taxa names
    first_uce_taxa <-  gsub(" .*","",first_uce[2:(numtaxa+1)])
    
    # Getting the first sequence line together
    seq_line <-  first_uce[2:(numtaxa+1)]
    for (j in 1:length(seq_line)) {
      seq_line[j] <- gsub(" ","",gsub(first_uce_taxa[j],"",seq_line[j]))
    }
    
    # Concatenating the rest of the sequences
    for (j in 1:length(seq_line)) {
      seq_line[j] <- gsub(" ","",paste(seq_line[j],paste(gsub(" ","",first_uce[seq((numtaxa+1+j),length(first_uce),numtaxa)]),collapse=""),collapse=""))
    }
    
  }
  
}