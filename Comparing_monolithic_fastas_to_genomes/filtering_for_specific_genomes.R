# Loading required libraries
library(tidyverse)

# Which genomes are required in order for a locus to be 
required_genomes <- c("Helmac","Quaclu")

# For completeness filtering step, how many taxa can a locus be missing from 
# and still incorporated in the final dataset?
allowed_taxa_drops <- 1

# Reading in non-paralogous loci
temp <- read_delim("/Users/alanaalexander/Dropbox/Hydrophilid_probe_design/blast_test/non_paralogous_uce_loci_blast_search.txt",delim="\t")

# Reading in original blast search to obtain the names of the blast genomes
ucegenome <- (read_delim("/Users/alanaalexander/Dropbox/Hydrophilid_probe_design/blast_test/blast_search_results.txt",delim=" ") %>% 
  select(genomename) %>% unique() %>% as.matrix())[,1]

# Creating a variable, full_spp, with all the genomes in a single column
paralogy_filtered <- as_tibble(cbind(temp,(
  unite(temp,full_spp,2:(length(ucegenome)+1)) %>% select(full_spp)
  )))

# Filtering for the first required genome
required_genome_filter <- paralogy_filtered %>% filter(grepl(required_genomes[1],full_spp))

# Filtering for any other required genomes
for (i in 2:length(required_genomes)) {
  required_genome_filter <- required_genome_filter %>% filter(grepl(required_genomes[i],full_spp))
}

print(paste("After filtering to retain",paste(required_genomes,collapse=","),dim(required_genome_filter)[1],"loci remain"))

# Filtering for completeness
required_genome_filter <- required_genome_filter %>% mutate(no_taxa=0)

for (i in ucegenome) {
  required_genome_filter <- required_genome_filter %>% mutate(no_taxa=ifelse((grepl(i,full_spp)),no_taxa+1,no_taxa))
}

required_genome_filter <- required_genome_filter %>% select(-full_spp)

required_genome_filter %>% group_by(no_taxa) %>% count()

# Filtering to completeness (allowing no more than allowed_taxa_drops of blast # genomes to not have the locus present)
required_genome_filter <- required_genome_filter %>% filter(no_taxa>(length(ucegenome)-allowed_taxa_drops))

write_delim(required_genome_filter,paste("retained_",paste(required_genomes,collapse="_"),"_completeness_filtered_uce_loci_atd_",allowed_taxa_drops,".txt",sep=""),col_names = TRUE,delim="\t")
print(paste("There are ",dim(required_genome_filter)[1], " UCE loci that were missing from no more than ",allowed_taxa_drops, " genomes and found in ",paste(required_genomes,collapse="_"),sep=""))
print(paste("This file has been written to retained_",paste(required_genomes,collapse="_"),"_completeness_filtered_uce_loci_atd_",allowed_taxa_drops,".txt",sep=""))

