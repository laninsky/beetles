# Comparing monolithic fastas to genomes
[comparing_monolithic_UCE_fastas](https://github.com/laninsky/comparing_monolithic_UCE_fastas) compares different monolithic fasta files containing all UCE loci for all taxa, identified using different base genomes in order to select the "best" base genome for probe design.

However, if your genomes are large, the genome alignment step required for comparing_monolithic_UCE_fastas becomes prohibitively long. These scripts implement an alternative where the existing list of loci in a single monolithic fasta file (rather than these having to be constructed for each genome) is blasted to each assembled genomes. If a given locus "hits" on multiple different regions on a given genome, they are marked as paralagous. However, if the loci make it through this first step, the next step is extracting each locus from the genome and then blasting it against our list of putative UCE loci (i.e. a reciprocal blast). If its best match is back to the locus that was originally used to identify it, than we can be pretty confident there is a one to one match between that UCE locus and that genome.

This process is repeated for all UCE loci and all genomes.

I had shared with you via dropbox a folder called Hydrophilid_probe_design that should have within:

(1) all the genomes in both fasta and 2bit format
(2) the monolithic.fasta file constructed during probe design
(3) the hydrophilid uce probes (I think there is a duplicate of this file in there with two different names)


# Need to make databases for each of our subject species
module load BLAST/2.9.0-gimkl-2018b
makeblastdb -in Rrat.faa -dbtype prot 
makeblastdb -in Rnorproteinseq.faa -dbtype prot 
makeblastdb -in guineapig_proteinseq.faa -dbtype prot 
makeblastdb -in goldenhamster_proteinseq.faa -dbtype prot 
makeblastdb -in ferret_proteinseq.faa -dbtype prot 
makeblastdb -in chinesehamster_proteinseq.faa -dbtype prot 

# Queried mm_queryproteinlist.faa against each of the
# subject genomes to find intial matches
# Called script Initial_blast.sh
# Submit by "sbatch Initial_blast.sh"

#!/bin/bash -e

#SBATCH -A uoo02820
#SBATCH -J Initial_blast
#SBATCH --ntasks 1
#SBATCH -c 1
#SBATCH -t 0:30:00
#SBATCH --mem=3G
#SBATCH -D /nesi/nobackup/uoo02820/manual_blast_code 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anna.clark@postgrad.otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread
#SBATCH --partition=large

module load BLAST/2.9.0-gimkl-2018b

blastp -query mm_queryproteinlist.faa -db Rrat.faa -evalue 1e-15 -outfmt 10 > Rrat_initial_blast.txt

blastp -query mm_queryproteinlist.faa -db Rnorproteinseq.faa -evalue 1e-15 -outfmt 10 > Rnorproteinseq_initial_blast.txt

blastp -query mm_queryproteinlist.faa -db guineapig_proteinseq.faa -evalue 1e-15 -outfmt 10 > guineapig_proteinseq_initial_blast.txt

blastp -query mm_queryproteinlist.faa -db goldenhamster_proteinseq.faa -evalue 1e-15 -outfmt 10 > goldenhamster_proteinseq_initial_blast.txt

blastp -query mm_queryproteinlist.faa -db ferret_proteinseq.faa -evalue 1e-15 -outfmt 10 > ferret_proteinseq_initial_blast.txt

blastp -query mm_queryproteinlist.faa -db chinesehamster_proteinseq.faa -evalue 1e-15 -outfmt 10 > chinesehamster_proteinseq_initial_blast.txt

# Summarizing the results of the previous blast searches and creating
# output *.faa files to then search back against mouse
# Called script initial_blast_filter.sh
# Submit by "sbatch initial_blast_filter.sh"

#!/bin/bash -e

#SBATCH -A uoo02820
#SBATCH -J Initial_blast_filter
#SBATCH --ntasks 1
#SBATCH -c 1
#SBATCH -t 0:30:00
#SBATCH --mem=3G
#SBATCH -D /nesi/nobackup/uoo02820/manual_blast_code 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anna.clark@postgrad.otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread
#SBATCH --partition=large

module load R/3.6.1-gimkl-2018b 

Rscript initial_blast_filter.R

# Now, need to carry out reciprocal blast
# First step is to database-ify our total mouse protein file
module load BLAST/2.9.0-gimkl-2018b
makeblastdb -in mmproteinrefseq.faa -dbtype prot 

# Then need to match our initial blast files back to this database
# Called submission script Reciprocal_blast.sh
# submitted it by sbatch Reciprocal_blast.sh

#!/bin/bash -e

#SBATCH -A uoo02820
#SBATCH -J Reciprocal_blast
#SBATCH --ntasks 1
#SBATCH -c 1
#SBATCH -t 0:30:00
#SBATCH --mem=3G
#SBATCH -D /nesi/nobackup/uoo02820/manual_blast_code 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anna.clark@postgrad.otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread
#SBATCH --partition=large

module load BLAST/2.9.0-gimkl-2018b

for i in *initial_blast.faa;
  do outputname=`echo $i | sed 's/_initial_blast.faa//g'`;
  blastp -query $i -db mmproteinrefseq.faa -evalue 1e-15 -outfmt 10 > $outputname.reciprocal_blast.txt;
done

# Our final-ish step will be to summarize reciprocal blast results between our
# species for the mouse proteins of interest
# Called submission script Reciprocal_summary.sh
# submitted it by sbatch Reciprocal_summary.sh

#!/bin/bash -e

#SBATCH -A uoo02820
#SBATCH -J Initial_blast_filter
#SBATCH --ntasks 1
#SBATCH -c 1
#SBATCH -t 0:30:00
#SBATCH --mem=3G
#SBATCH -D /nesi/nobackup/uoo02820/manual_blast_code 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anna.clark@postgrad.otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread
#SBATCH --partition=large

module load R/3.6.1-gimkl-2018b 

Rscript reciprocal_blast_summary.R



