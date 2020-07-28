# Comparing monolithic fastas to genomes
[comparing_monolithic_UCE_fastas](https://github.com/laninsky/comparing_monolithic_UCE_fastas) compares different monolithic fasta files containing all UCE loci for all taxa, identified using different base genomes in order to select the "best" base genome for probe design.  

However, if your genomes are large, the genome alignment step required for comparing_monolithic_UCE_fastas becomes prohibitively long. These scripts implement an alternative where the existing list of loci in a single monolithic fasta file (rather than these having to be constructed for each genome) is blasted to each assembled genome (see code further down in this file).

After this step, initial_blast_filter.R is used to determine if a given locus/monolithic taxon combo has "good" hits to multiple different regions on a given genome. If so, they are marked as paralagous. If different versions of the uce-locus (i.e. deriving from different monolithic taxa) have good hits to different regions in a given genome, they are also marked as paralogous.

To be considered a good match, the matches have to be more than 120 bp (this cut-off is used as the normal length of the probe sequence for uce-loci is 120 bp. It can be modified by finding and replacing "120" in initial_blast_filter.R) and/or have an e-value of 0. If there is only one match more than 120 bp and/or with an e-value of 0, the locus is recorded as "good". If there are no matches more than 120 bp and/or with an e-value, that uce_taxon is recorded as "NO_STRONG_MATCH" for that given blast'ed genome.  If there are multiple "good" matches the locus is recorded as paralogous as described above.

This process is repeated for all UCE loci and all genomes.  

### What you need:
* Your \*.monolithic.fasta file  
* Your genomes in fasta format in the same location as the \*.monolithic.fasta.file  
* BLAST to be in your path (e.g. export PATH=/Users/alanaalexander/Dropbox/Hydrophilid_probe_design/blast_test/bin/ncbi-blast-2.10.0+/bin:$PATH)
* R in your path
* The R package tidyverse to be installed

### Options to be set before running
First, you need to set the similarity threshold you want to use for your BLAST searches:  
`export BLASTSIM=99`

You also need to modify the 6th line of initial_blast_filter.R to specify the number of blast'ed
genomes a uce-locus can be missing in but still be incorporated in the dataset. The default setting is 2 (e.g. a locus has to be not found in a blast search for more than two of your total genomes before being excluded from the dataset).

### In your folder with your files:
Make blast databases for each fasta file (i.e. the monolithic file, and your genomes)
```
for i in *.fasta;
  do makeblastdb -in $i -dbtype nucl;
done
```
Specify the monolithic file versus the genome files
```
monolithic_file=`ls *.fasta | grep 'monolithic'`
genome_files=`ls *.fasta | grep -v 'monolithic'`
```
Run the initial blast search
```
for i in $genome_files;
  do blastn -db $i -query $monolithic_file -perc_identity $BLASTSIM -outfmt '10 sseqid qseqid length qlen sstart send evalue' > $i.initialblast.txt;
done
```
Conduct the comparison across genomes
```
Rscript initial_blast_filter.R
```

### Output files:
From the initial blast searches  
* blast database files for each of the genomes and monolithic fasta files
* \*.initialblast.txt files for each genome, giving the blast matches between the genome and the monolithic fasta file
From initial_blast_filter.R, a series of tab-delimited files:
* blast_search_results.txt: a summary of the blast results within genomes, including whether loci/monolithic taxon combos were considered "PARALOGOUS" or "NO_STRONG_MATCH" in each genome (see above for descriptions of each of these categories). Columns include 'genomename' (the genome that was the subject of the blast search), blast output, and the uce_locus and uce_genome (taxa from the monolithic fasta that the uce locus was characterized in).
* full_uce_loci_blast_search.txt: a summary after comparing across ghe blast genomes. There is one row per uce_locus (column uce_loci), and then a column for each of the monolithic taxa that these loci were characterized in. The contents of these columns are the blast genomes where a match was found for that uce_locus/monolithic taxon combination, separated by comma if multiple genomes had matches to that uce_locus/monolithic taxon combination. If any of these matches were in the "PARALOGOUS" or "NO_STRONG_MATCH" categories, this is also noted. A final set of columns is prefaced by "length_", corresponding to each of these monolithic taxa. The contents of these columns are the average length of the matches across the genomes that had good blast mathches to that uce_loci/monolithic taxon combination.
* non_paralogous_uce_loci_blast_search.txt: The subset of loci from full_uce_loci_blast_search.txt that made it through the paralogy filter.
* paralogy_filtered_uce_loci.txt: The subset of loci from full_uce_loci_blast_search.txt that didn't make it through the paralogy filter.
* completeness_filtered_uce_loci_atd_X.txt, where X = the number you specify for allowed_taxa_drops in initial_blast_filter.R line 6: The subset of loci from non_paralogous_uce_loci_blast_search.txt that were missing "good" blast matches to no more than X genomes.

If initial files in the R pipeline have already been created in the directory this analysis is conducted in, they will not be recreated but read in from the directory instead.

### Futher filtering
filtering_for_specific_genomes.R allows you to specify on Line 5 any genomes that are required to be represented in the loci that make it through your filters. You can then also implement allowed_taxa_drops on Line 9 if you want to additionally filter for completeness (default 2). You also need to modify the location of your non_paralogous_uce_locu_blast_search.txt file on Line 12, and the location of your blast_search_results.txt file on Line 15
