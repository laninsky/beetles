# Load libraries that we need
library(tidyverse)

# For filtering step, how many taxa can a locus be missing from and still
# incorporated in the final dataset?
allowed_taxa_drops <- 1

# If the pipeline has been run before and blast_search_results.txt exists in the directory
# then reading it in
if ("blast_search_results.txt" %in% list.files()) {
  output <- read_table2("blast_search_results.txt")
} else {
  # Creating the output file that summarizes the blast searches  
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

  # Writing out the summary of the blast search results
  write_delim(output,"blast_search_results.txt",col_names = TRUE)
  print("The blast search results have been written to blast_search_results.txt")
}

# Getting the total list of uce loci
ucelocilist <- unique(output$uce_locus)
# Getting the total list of taxa the monolithic file had loci characterized over
ucegenome <- unique(output$uce_genome)
# Getting the list of the genomes that were blasted
genomes <-  unique(output$genomename)

# If the pipeline has been run before and full_uce_loci_blast_search.txt exists in the directory
# then reading it in
if ("full_uce_loci_blast_search.txt" %in% list.files()) {
  evaluating_loci <- read_delim("full_uce_loci_blast_search.txt",delim="\t")
} else {
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
  write_delim(evaluating_loci,"full_uce_loci_blast_search.txt",col_names = TRUE,delim="\t")
  print(paste("There are ",dim(evaluating_loci)[1], " total UCE loci developed across the following taxa:",sep=""))
  print(ucegenome)
  print("That were found in blast searches to the following genomes:")
  print(genomes)
  print("The full list of these loci has been written to full_uce_loci_blast_search.txt")
}

# If the pipeline has been run before and non_paralogous_uce_loci_blast_search.txt and
# paralogy_filtered_uce_loci.txt exist in the directory then reading them in
if ("non_paralogous_uce_loci_blast_search.txt" %in% list.files() & "paralogy_filtered_uce_loci.txt" %in% list.files()) {
  paralogy_filtered <- read_delim("non_paralogous_uce_loci_blast_search.txt", delim="\t")
  loci_filtered_due_to_paralogy <- read_delim("paralogy_filtered_uce_loci.txt",delim="\t")
} else {

  # Filtering out loci that show paralogy in any of the single blast searches
  paralogy_filtered <- evaluating_loci
  loci_filtered_due_to_paralogy <- NULL

  for (i in ucegenome) {
    paralogy_filtered <- paralogy_filtered %>% filter(!grepl("PARALOGOUS",!!as.name(i)))
    loci_filtered_due_to_paralogy <- rbind(loci_filtered_due_to_paralogy,(evaluating_loci %>% filter(grepl("PARALOGOUS",!!as.name(i)))))
  }
  
  surviving_loci_count <- dim(paralogy_filtered)[1]
  paralogy_within_count <- dim(loci_filtered_due_to_paralogy)[1]
  
  print(paste(paralogy_within_count,"uce loci were removed from the dataset because of paralogous matches from a single uce locus/monolithic taxon combination to a genome."))
  print(paste(surviving_loci_count,"uce loci made it through this filter"))
  print(paste("Now filtering for uce loci that are paralogous based on blast searches when using different monolithic taxon versions of a given uce locus"))
  
  # Filtering out loci that are paralogous across the different monolithic taxa versions of each uce locus
  for (i in 1:length(genomes)) {
    for (j in 1:length(ucelocilist)) {
      temp <- output %>% filter(genomename==genomes[i]) %>% filter(uce_locus==ucelocilist[j]) %>% filter(sseqid!="NO_STRONG_MATCH")
      # If more than one monolithic taxa had a blast match to that genome for that uce locus...
      if (dim(temp)[1]>1) {
        # If the different uce_taxon for a given uce matches to a different scaffold in the genome
        if (all(temp$sseqid!=temp$sseqid[1])) {
          print(paste("Removing",ucelocilist[j],"as it is paralogous across blasts of different monolithic taxa versions for",genomes[i]))
          # Adding this loci to the loci_filtered due to paralogy dataset
          loci_filtered_due_to_paralogy <- rbind(loci_filtered_due_to_paralogy,(paralogy_filtered %>% filter(uce_loci==ucelocilist[j])))
          # And removing this uce locus from the paralogy filtered dataset
          paralogy_filtered <- paralogy_filtered %>% filter(uce_loci!=ucelocilist[j])
        } else {
          # Double checking that the location on the scaffold of the different matches is consistent
          for (k in 2:dim(temp)[1]) {
            # If the ranges for any do not overlap, writing these loci out, and then breaking out of the loop
            if (!((max(temp$sstart[1],temp$send[1])>min(temp$sstart[k],temp$send[k])) & ((min(temp$sstart[1],temp$send[1])<max(temp$sstart[k],temp$send[k]))))) {
              print(paste("Removing",ucelocilist[j],"as it is paralogous across blasts of different monolithic taxa versions for",genomes[i]))
              # Adding this loci to the loci_filtered due to paralogy dataset
              loci_filtered_due_to_paralogy <- rbind(loci_filtered_due_to_paralogy,(paralogy_filtered %>% filter(uce_loci==ucelocilist[j])))
              # And removing this uce locus from the paralogy filtered dataset
              paralogy_filtered <- paralogy_filtered %>% filter(uce_loci!=ucelocilist[j])
              break
            }
          }
        }
      } 
    }
  }
  
  print(paste("In total, an additional",(dim(loci_filtered_due_to_paralogy)[1]-paralogy_within_count), "loci were removed"))
  
  # Writing out the summary of all the paralogy filter
  write_delim(paralogy_filtered,"non_paralogous_uce_loci_blast_search.txt",col_names = TRUE, delim="\t")
  write_delim(loci_filtered_due_to_paralogy,"paralogy_filtered_uce_loci.txt",col_names = TRUE, delim="\t")
  print(paste("There are ",dim(loci_filtered_due_to_paralogy)[1], " UCE loci that were filtered due to paralogy",sep=""))
  print(paste("There are ",dim(paralogy_filtered)[1], " UCE loci that made it through the paralogy filter",sep=""))
  print("The paralogous loci have been written to paralogy_filtered_uce_loci.txt")
  print("The nonparalogous loci have been written to non_paralogous_uce_loci_blast_search.txt")
}
  
  
completeness_filter <- as_tibble(cbind(paralogy_filtered,(paralogy_filtered %>% unite(full_spp,2:length(ucegenome)) %>% select(full_spp)))) %>% mutate(no_taxa=0)

for (i in ucegenome) {
  completeness_filter <- completeness_filter %>% mutate(no_taxa=ifelse((grepl(i,full_spp)),no_taxa+1,no_taxa))
}

completeness_filter <- completeness_filter %>% select(-full_spp)

# If the pipeline has been run before and blast_search_results.txt exists in the directory
# then reading it in
if (paste("completeness_filtered_uce_loci_atd_",allowed_taxa_drops,".txt",sep="") %in% list.files()) {
  output <- read_table2(paste("completeness_filtered_uce_loci_atd_",allowed_taxa_drops,".txt",sep=""))
} else {
  # Filtering to completeness (allowing no more than allowed_taxa_drops of blast genomes
  # to not have the locus present)
  completeness_filter <- completeness_filter %>% filter(no_taxa>(length(genomes)-allowed_taxa_drops))
  write_delim(completeness_filter,paste("completeness_filtered_uce_loci_atd_",allowed_taxa_drops,".txt",sep=""),col_names = TRUE,delim="\t")
  print(paste("There are ",dim(completeness_filter)[1], " UCE loci that were missing from no more than ",allowed_taxa_drops, " genomes that were blast searched",sep=""))
  print(paste("These loci have been written to completeness_filtered_uce_loci_atd_",allowed_taxa_drops,".txt",sep=""))
}

