# Comparing monolithic fastas to genomes
[comparing_monolithic_UCE_fastas](https://github.com/laninsky/comparing_monolithic_UCE_fastas) compares different monolithic fasta files containing all UCE loci for all taxa, identified using different base genomes in order to select the "best" base genome for probe design.  

However, if your genomes are large, the genome alignment step required for comparing_monolithic_UCE_fastas becomes prohibitively long. These scripts implement an alternative where the existing list of loci in a single monolithic fasta file (rather than these having to be constructed for each genome) is blasted to each assembled genomes. If a given locus "hits" on multiple different regions on a given genome, they are marked as paralagous. However, if the loci make it through this first step, the next step is extracting each locus from the genome and then blasting it against our list of putative UCE loci (i.e. a reciprocal blast). If its best match is back to the locus that was originally used to identify it, than we can be pretty confident there is a one to one match between that UCE locus and that genome.  

This process is repeated for all UCE loci and all genomes.  

### What you need:
* Your \*.monolithic.fasta file  
* Your genomes in fasta format in the same location as the \*.monolithic.fasta.file  
* BLAST to be in your path (e.g. export PATH=/Users/alanaalexander/Dropbox/Hydrophilid_probe_design/blast_test/bin/ncbi-blast-2.10.0+/bin:$PATH)  

### Options to be set before running
First, you need to set the similarity threshold you want to use for your BLAST searches:  
`export BLASTSIM=99`

### In your folder with your files:
```
# Make blast databases for each fasta file
for i in *.fasta;
  do makeblastdb -in $i -dbtype nucl;
done

# Specifying the monolithic file versus the genome files
monolithic_file=`ls *.fasta | grep 'monolithic'`
genome_files=`ls *.fasta | grep -v 'monolithic'`

# Running the initial blast search
for i in $genome_files;
  do blastn -db $i -query $monolithic_file -perc_identity $BLASTSIM -outfmt '10 sseqid qseqid length qlen sstart send evalue' > $i.initialblast.txt;
done

# Modify the 6th line of initial_blast_filter.R to specify the number of blast'ed
# genomes a uce-locus can be missing in but still be incorporated in the dataset
Rscript initial_blast_filter.R



