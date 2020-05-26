# Load libraries that we need
library(tidyverse)

# Get list of initial_blast results
initialblasts <- list.files(pattern="initial_blast.txt")

# For each of our initial_blast results
for (i in initialblasts) {
  print(paste("Up to file ",i,sep=""))
  
  # Getting corresponding faa file name
  tempfaaname <- gsub("_initial_blast.txt",".faa",i)
  
  # Reading in our blast results
  tempblast <- read_delim(i,delim=",",col_names = FALSE)
  
  # Take best hit(s) for each mouse scaffold
  besthits <- tempblast %>% group_by(X1) %>% filter(X11==min(X11)) %>% filter(X12==max(X12))
  
  # Our subject faa files have isoforms, so we need
  # to create a new variable where all the isoforms
  # have the same name
  besthits <- besthits %>% mutate(gene_name=gsub("-.*","",X2))
  
  # Summarize results where multiple hits found
  # (these will include different isoforms)
  print("Summarizing multiple hits")
  
  multiplehitloci <- besthits %>% group_by(X1) %>% summarise(count=n()) %>% filter(count>1) %>% select(X1) %>% as.matrix()
  
  isoforms <- NULL
  multiloci <- NULL
  # Sorting through this to see if isoforms or not
  print("Distinguishing between multiple genes and multiple isoforms per gene")
  
  for (j in multiplehitloci) {
    temp <- besthits %>% filter(X1==j)
    if(dim(temp %>% group_by(gene_name) %>% summarise(count=n()))[1]==1) {
      unique(temp$gene_name)
     isoforms <- c(isoforms,unique(temp$gene_name))
   } else {
     multiloci <- c(multiloci,unique(temp$gene_name))
   }
  }
  
  # Write out the multiple hit loci to a file
  print("Writing out multiple hit loci to a multiple_hit_loci.txt")
  if (!is.null(multiloci)) {
    multilociout <- besthits %>% filter(X1 %in% multiplehitloci) %>% mutate(genome=i)
    write_delim(multilociout,"multiple_hit_loci.txt",append=TRUE)
  
  # Now back to our non-multiple-hit loci
    print("Summarizing our non-multiple best hit loci")
    nonmultiloci <- besthits %>% filter(!(X1 %in% multiplehitloci)) %>% mutate(genome=i)
  } else {
    nonmultiloci <- besthits
  }
  
  # Now we need to deal with any isoforms
  print("Writing out isoforms to isoform_record.txt, and taking one representative isoform per loci")
  if (!is.null(isoforms)) {
    nonmultiloci <- nonmultiloci %>% arrange(X1,gene_name,desc(X4),desc(X3))
    for (j in isoforms) {
      temp <- nonmultiloci %>% filter(gene_name==j) %>% mutate(genome=i)
      write_delim(temp,"isoform_record.txt",append=TRUE)
      toremove <- temp$X2[which(!(temp$X2 %in% temp$X2[1]))]
      nonmultiloci <- nonmultiloci %>% filter(!(X2 %in% toremove))      
      
    }
  }
  
  # First, need to check that multiple subject loci
  # do not match to different query loci
  print("Writing out subject loci that match to different query loci: duplicated_in_query.txt")
  duplicated_in_query <- nonmultiloci %>% group_by(X2) %>% summarise(count=n()) %>% filter(count>1) %>% select(X2) %>% as.matrix()
  
  duplloci <- nonmultiloci %>% filter(X2 %in% duplicated_in_query) %>% mutate(genome=i)

  write_delim(duplloci,"duplicated_in_query.txt",append=TRUE)
  
  # So filter any nonmultiloci that are actually duplloci
  nonmultiloci <- nonmultiloci %>% filter(!(X2 %in% duplicated_in_query)) %>% mutate(genome=i) 
  
  print("Writing out 'good loci': non_multiple_hit_loci.txt")
  write_delim(nonmultiloci,"non_multiple_hit_loci.txt",append=TRUE)
  
  # Filter faa file based on best hit for reciprocal blast
  # First, need to read in *.faa file
  tempfaa <- read_delim(tempfaaname,col_names = FALSE,delim="@")
  
  # Getting our subject scaffold names and our newnames for these
  subj_scaffolds <- paste(">",nonmultiloci$X2,sep="")
  newnames <- paste(paste(">",nonmultiloci$X2,sep=""),nonmultiloci$X1,sep="___")
  
  # Where are all the headerlines in the *.faa file?
  headerlines <- which(grepl(">",tempfaa$X1))
  
  # Stripping out the descriptor after the sequence names in the faa files
  # to match the blast results
  tempfaa[headerlines,1] <- gsub(" .*","",as.matrix(tempfaa[headerlines,1]))
  
  # Where in faa do our subj_scaffold names occur?
  scaffheaderlines <- which(tempfaa$X1 %in% subj_scaffolds)
  
  # What position do the sequences end on
  endpos <- headerlines[which(headerlines %in% scaffheaderlines)+1]-1
  
  # Name for our outfile
  tempoutfile <- gsub(".txt",".faa",i)
  
  # Generating our outfile
  temp <- NULL
  for (j in 1:length(scaffheaderlines)) {
    temp <- bind_rows(temp,tempfaa[scaffheaderlines[j]:endpos[j],1])
  }
  
  # Replacing exisiting sequence name with preferred name
  namestoreplace <- which(grepl(">",temp$X1))
  
  for (j in namestoreplace) {
    temp[j,1] <- newnames[which(grepl(temp$X1[j],newnames))]
  }
  
  print("Writing out fasta file with sequences for reverse blast")
  write_delim(temp,tempoutfile,append=TRUE)
  
}