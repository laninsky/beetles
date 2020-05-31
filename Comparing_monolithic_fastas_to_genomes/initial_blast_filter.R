# Load libraries that we need
library(tidyverse)

# For filtering step, how many taxa can a locus be missing from and still
# incorporated in the final dataset?
allowed_taxa_drops <- 2

# Get list of initial_blast results
initialblasts <- list.files(pattern="initialblast.txt")

output <- NULL

# For each of our initial_blast results
for (i in initialblasts) {
  print(paste("Up to file ",i,sep=""))
  
  # Getting corresponding genome name
  genomename <- gsub(".fasta.initialblast.txt","",i)
  
  # Reading in our blast results for that genome
  tempblast <- read_delim(i,delim=",",col_names = FALSE)
  names(tempblast) <- c("sseqid","qseqid","length","qlen","sstart","send","evalue")
  
  # Getting the unique qseqids (locus_taxon combinations) in the initialblast.txt file
  uniqueqseqid <- unique(tempblast$qseqid)
  
  # Tempoutput variable for filtering within a locus_taxon
  # for too short matches and paralagous loci
  tempoutput <- NULL
  
  # Stepping through these unique qseqids and filtering
  for (j in uniqueqseqid) {
    # Filtering to results for each locus_taxon combination
    temp <- tempblast %>% filter(qseqid==j) 
    # If only one blast result returned
    if (dim(temp)[1]==1) {
      # If the length is greater than 120 bp (the average length of probes)
      if (temp$length > 120) {
        # Then record this match as "good"
        tempoutput <- rbind(tempoutput,temp)
      }
    } else {
      # If there is more than one match, filtering for evalue==0 and length > 120 bp
      temp <- temp %>% filter(evalue==0 & length > 120)
      # If there are one or more matches fulfilling this criteria
      if (dim(temp)[1]>=1) {
        # If there is just one match that fulfils this criteria
        if (dim(temp)[1]==1) {
          # Then record this match as "good"
          tempoutput <- rbind(tempoutput,temp)
        } else {
          # If there are multiple "good" matches, recording this uce_locus as paralogous
          temp <- temp[1,]
          temp$sseqid <- "PARALOGOUS"
          tempoutput <- rbind(tempoutput,temp)
        }  
      } else {
        # If there weren't any matches with an evalue==0 and a length > 120
        # then just filtering for evalue
        temp <- tempblast %>% filter(qseqid==j) %>% filter(evalue==0)
        # If only one match retained
        if (dim(temp)[1]==1) {
          # Then writing this match out as good
          tempoutput <- rbind(tempoutput,temp)
        } else {
          # If 0 or more than one match had an evalue of 0, filtering for length instead
          temp <- tempblast %>% filter(qseqid==j) %>% filter(length > 120)
          # If one or more match recorded
          if (dim(temp)[1]>=1) {
            # If only one match fulfils this criteria
            if (dim(temp)[1]==1) {
              # Then this locus will be written out as "good"
              tempoutput <- rbind(tempoutput,temp)
            } else {
              # If multiple loci returned, then recording them as paralogous
              temp <- temp[1,]
              temp$sseqid <- "PARALOGOUS"
              tempoutput <- rbind(tempoutput,temp)
            }  
          } else {
            # If one or more matches aren't returned that are greater than 120 bp
            # Then filtering e-values not being zero and length being less than 120 bp
            temp <- tempblast %>% filter(qseqid==j)  %>% filter(length < 120 | evalue!=0)
            # Writing this locus out as not having a strong match
            temp <- temp[1,]
            temp$sseqid <- "NO_STRONG_MATCH"
            tempoutput <- rbind(tempoutput,temp)
          }  
        }
      }  
    }
  }
  # Recording this locus with the genome it was blasted in
  tempoutput <- as_tibble(cbind(genomename,tempoutput))
  # Recording this output in the "final" output
  output <- rbind(output,tempoutput)
}  

# Splitting the uce_locus/taxon it was characterized in for the monolithic file
# into separate columns
output <- output %>% mutate(uce_locus=gsub("_.*","",qseqid)) %>% mutate(uce_genome=gsub(".*_","",qseqid))

# Getting the total list of uce loci
ucelocilist <- unique(output$uce_locus)
# Getting the total list of taxa the monolithic file had loci characterized over
ucegenome <- unique(output$uce_genome)
# Getting the list of the genomes that were blasted
genomes <-  unique(output$genomename)

# Creating an output matrix to hold our comparisons across genomes
evaluating_loci <- matrix(NA,ncol=length(ucegenome)*2,nrow=length(ucelocilist))

# For each genome the monolithic file had loci characterized over
for (i in 1:length(ucegenome)) {
  # For each loci
  for (j in 1:length(ucelocilist)) {
    # Filtering for whether any blast results were paralogous for that combination of
    # uce locus and taxon it was characterized in in the monolithic file
    temp <- output %>% filter(uce_locus==ucelocilist[j]) %>% filter(uce_genome==ucegenome[i]) %>% filter(sseqid=="PARALOGOUS")
    # If there is more than one paralogous result
    if (dim(temp)[1]>0) {
      # Then recording for this uce locus/monolithic taxon the blast genomes that
      # were found to have paralogous loci
      evaluating_loci[j,i] <- paste("PARALOGOUS in",paste(as.matrix(temp[,1]),collapse=","))
    } else {
      # Otherwise filtering to "actual" results for that combination of
      # uce locus and genome it was characterized in in the monolithic file
      temp <- output %>% filter(uce_locus==ucelocilist[j]) %>% filter(uce_genome==ucegenome[i]) %>% filter(sseqid!="PARALOGOUS") %>% filter(sseqid!="NO_STRONG_MATCH")
      # If there is at least one result for this
      if (dim(temp)[1]>0) {
        # Recording the blast genomes this locus/monolithic taxon was found in
        evaluating_loci[j,i] <- paste(as.matrix(temp[,1]),collapse=",")
        evaluating_loci[j,(length(ucegenome)+i)] <- mean(temp$length,na.rm=TRUE)
      }
    }
  }
}

# Attaching the uce names to the result matrix
evaluating_loci <- as_tibble(cbind(ucelocilist,evaluating_loci))
# Giving names to the columns (which refer to the taxon that loci in the monolithic
# file were characterized in. Contents of these columns refer to the blast genomes
# these loci were found in
names(evaluating_loci) <- c("uce_loci",ucegenome,paste("length_",ucegenome,sep=""))

# Writing out the summary of all the uce-loci in the blast search
write_delim(evaluating_loci,"full_uce_loci_blast_search.txt",col_names = TRUE)
print(paste("There are ",dim(evaluating_loci)[1], " total UCE loci developed across the following taxa:",sep=""))
print(ucegenome)
print("That were found in blast searches to the following genomes:")
print(genomes)
print("The full list of these loci has been written to full_uce_loci_blast_search.txt")

# Filtering out loci that show paralogy in any of the blast searches
paralogy_filtered <- evaluating_loci
loci_filtered_due_to_paralogy <- NULL

for (i in ucegenome) {
  paralogy_filtered <- paralogy_filtered %>% filter(!grepl("PARALOGOUS",!!as.name(i)))
  loci_filtered_due_to_paralogy <- rbind(loci_filtered_due_to_paralogy,(evaluating_loci %>% filter(grepl("PARALOGOUS",!!as.name(i)))))
}

# Writing out the summary of all the paralogy filter
write_delim(paralogy_filtered,"non_paralogous_uce_loci_blast_search.txt",col_names = TRUE)
write_delim(loci_filtered_due_to_paralogy,"paralogy_filtered_uce_loci.txt",col_names = TRUE)
print(paste("There are ",dim(loci_filtered_due_to_paralogy)[1], " UCE loci that were filtered due to paralogy",sep=""))
print(paste("There are ",dim(paralogy_filtered)[1], " UCE loci that made it through the paralogy filter",sep=""))
print("The paralogous loci have been written to paralogy_filtered_uce_loci.txt")
print("The nonparalogous loci have been written to non_paralogous_uce_loci_blast_search.txt")

completeness_filter <- as_tibble(cbind(paralogy_filtered,(paralogy_filtered %>% unite(full_spp,2:length(ucegenome)) %>% select(full_spp)))) %>% mutate(no_taxa=0)

for (i in ucegenome) {
  completeness_filter <- completeness_filter %>% mutate(no_taxa=ifelse((grepl(i,full_spp)),no_taxa+1,no_taxa))
}

# Filtering to completeness (allowing no more than allowed_taxa_drops of blast genomes
# to not have the locus present)
completeness_filter <- completeness_filter %>% filter(no_taxa>(length(genomes)-allowed_taxa_drops))
write_delim(completeness_filter,"completeness_filtered_uce_loci.txt",col_names = TRUE)
print(paste("There are ",dim(completeness_filter)[1], " UCE loci that were missing from no more than ",allowed_taxa_drops, " genomes that were blast searched",sep=""))


