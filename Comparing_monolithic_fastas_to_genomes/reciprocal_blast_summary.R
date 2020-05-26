# Load libraries that we need
library(tidyverse)

# Get list of initial_blast results
recipblasts <- list.files(pattern=".reciprocal_blast.txt")

# Create our cross_species variable that will hold out comparison
# across species
cross_species <- NULL

# For each of our initial_blast results
for (i in recipblasts) {
  print(paste("Up to file ",i,sep=""))
  
  # Reading in our blast results
  tempblast <- read_delim(i,delim=",",col_names = FALSE)
  
  # Take best hit(s) for each query scaffold
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
  
  print("Writing out multiple hit loci to a reciprocal_multiple_hit_loci.txt")
  
  if (!is.null(multiloci)) {
    multilociout <- besthits %>% filter(X1 %in% multiplehitloci) %>% mutate(genome=i)
    write_delim(multilociout,"reciprocal_multiple_hit_loci.txt",append=TRUE)
    
    # Now back to our non-multiple-hit loci
    print("Summarizing our non-multiple best hit loci")
    nonmultiloci <- besthits %>% filter(!(X1 %in% multiplehitloci)) %>% mutate(genome=i)
  } else {
    nonmultiloci <- besthits
  }
  
  # Now we need to deal with any isoforms
  print("Writing out isoforms to reciprocal_isoform_record.txt, and taking one representative isoform per loci")
  if (!is.null(isoforms)) {
    nonmultiloci <- nonmultiloci %>% arrange(X1,gene_name,desc(X4),desc(X3))
    for (j in isoforms) {
      temp <- nonmultiloci %>% filter(gene_name==j) %>% mutate(genome=i)
      write_delim(temp,"reciprocal_isoform_record.txt",append=TRUE)
      toremove <- temp$X2[which(!(temp$X2 %in% temp$X2[1]))]
      nonmultiloci <- nonmultiloci %>% filter(!(X2 %in% toremove))
    }
  }
  
  # First, need to check that multiple subject loci
  # do not match to different query loci
  print("Writing out subject loci that match to different query loci: duplicated_in_query.txt")
  duplicated_in_query <- nonmultiloci %>% group_by(X2) %>% summarise(count=n()) %>% filter(count>1) %>% select(X2) %>% as.matrix()
  
  if (length(duplicated_in_query)>0) {
    duplloci <- nonmultiloci %>% filter(X2 %in% duplicated_in_query) %>% mutate(genome=i)  
    write_delim(duplloci,"reciprocal_duplicated_in_query.txt",append=TRUE)
  }
  
  # So filter any nonmultiloci that are actually duplloci
  nonmultiloci <- nonmultiloci %>% filter(!(X2 %in% duplicated_in_query))  
  
  # Now, to see whether the loci do match back to the original match
  # Need to pull out mouse scaffold name from X1
  different_gene <- nonmultiloci %>% filter(gsub(".*___","",X1)!=X2) %>% select(X1) %>% as.matrix()
  
  # Write out these failures to recapture the original protein sequence,
  # despite only returning one best hit
  
  if (length(different_gene)>0) {
    recap_failure <- nonmultiloci %>% filter(X1 %in% different_gene) %>% mutate(genome=i)
    print("Writing out 'recapture failures': reciprocal_recap_failures.txt")
    write_delim(recap_failure,"reciprocal_recap_failures.txt",append=TRUE)
    
    # So filter any nonmultiloci that are actually recap failures
    nonmultiloci <- nonmultiloci %>% filter(!(X1 %in% different_gene)) %>% mutate(genome=i)
  }
  
  print("Writing out 'good loci': reciprocal_non_multiple_hit_loci.txt")
  write_delim(nonmultiloci,"reciprocal_non_multiple_hit_loci.txt",append=TRUE)
  
  # Saving for comparison across species
  cross_species <- bind_rows(cross_species,(nonmultiloci %>% mutate(genome=i)))
  
}  
  
# Writing out summaries of the data
perc_summary <- spread(cross_species[,c(-1,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13)],genome,X3)
write_delim(perc_summary,"percentage_match_summary.txt")
length_summary <- spread(cross_species[,c(-1,-3,-5,-6,-7,-8,-9,-10,-11,-12,-13)],genome,X4)
write_delim(length_summary ,"length_summary.txt")

  