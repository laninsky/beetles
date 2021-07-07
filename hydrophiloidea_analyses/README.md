### Installing phyluce
Load the Miniconda module:
```
module load Miniconda3/4.9.2
```
Obtain and install phyluce:
```
wget https://raw.githubusercontent.com/faircloth-lab/phyluce/v1.7.1/distrib/phyluce-1.7.1-py36-Linux-conda.yml
conda env create --file phyluce-1.7.1-py36-Linux-conda.yml --prefix /nesi/nobackup/uoo00105/beetles/conda
```
Activate conda environment:
```
conda activate /nesi/nobackup/uoo00105/beetles/conda
```
After finishing up:
```
conda deactivate
```

### Uploading data (50% complete matrix)
```
scp -r /Users/alanaalexander/Dropbox/beetles/hydrophilids_me/50perc_copied_from_harddrive/abyss_50perc_nexus mahuika:/nesi/nobackup/uoo00105/beetles
```

### Running Phyluce 
When logging in following the previous installation step:
```
module load Miniconda3/4.9.2
conda activate /nesi/nobackup/uoo00105/beetles/conda
```
Obtaining the complete list of taxa present across these files
```
for i in *.nexus;
  do tail -n+6 $i | head -n-2 | cut -d " " -f 1 >> taxa_list.txt;
done

sort taxa_list.txt | uniq | wc -l

mv taxa_list.txt ../
```
Create some of the necessary directories:
```
mkdir logs
```

Obtaining the 70% complete matrix from the 50% matrix
```
#!/bin/bash -e

#SBATCH -A uoo00105 
#SBATCH -J 70_perc_matrix
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 15:00
#SBATCH --mem=20G
#SBATCH -D /nesi/nobackup/uoo00105/beetles
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread
#SBATCH --qos=debug

module load Miniconda3/4.9.2
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate /nesi/nobackup/uoo00105/beetles/conda
phyluce_align_get_only_loci_with_min_taxa --alignments abyss_50perc_nexus --taxa 63 --percent 0.7 --output 70perc_nexus --cores 12 --log-path logs
```

Obtaining the concatenated files for raxml and exabayes
```
#!/bin/bash -e

#SBATCH -A uoo00105 
#SBATCH -J 50perc_raxml
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 15:00
#SBATCH --mem=20G
#SBATCH -D /nesi/nobackup/uoo00105/beetles
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread
#SBATCH --qos=debug

module load Miniconda3/4.9.2
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate /nesi/nobackup/uoo00105/beetles/conda
phyluce_align_concatenate_alignments --alignments /nesi/nobackup/uoo00105/beetles/abyss_50perc_nexus --output /nesi/nobackup/uoo00105/beetles/50perc_raxml --phylip
```
```
#!/bin/bash -e

#SBATCH -A uoo00105 
#SBATCH -J 50perc_raxml
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 15:00
#SBATCH --mem=20G
#SBATCH -D /nesi/nobackup/uoo00105/beetles
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread
#SBATCH --qos=debug

module load Miniconda3/4.9.2
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate /nesi/nobackup/uoo00105/beetles/conda
phyluce_align_concatenate_alignments --alignments /nesi/nobackup/uoo00105/beetles/70perc_nexus --output /nesi/nobackup/uoo00105/beetles/70perc_raxml --phylip
```
### Making trees
Running raxml on the 70% files (50% had previously been run)
```
#!/bin/bash -e

#SBATCH -A uoo00105 
#SBATCH -J 70perc_run1
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 1:00:00
#SBATCH --mem=20G
#SBATCH -D /nesi/nobackup/uoo00105/beetles/70perc_raxml
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread

module load RAxML/8.2.12-gimkl-2020a
raxmlHPC-PTHREADS-SSE3 -s 70perc_raxml.phylip -n run1 -m GTRCATI -f a -N 100 -x $RANDOM -p $RANDOM -T 12
```
```
#!/bin/bash -e

#SBATCH -A uoo00105 
#SBATCH -J 70perc_run2
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 1:00:00
#SBATCH --mem=20G
#SBATCH -D /nesi/nobackup/uoo00105/beetles/70perc_raxml
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread

module load RAxML/8.2.12-gimkl-2020a
raxmlHPC-PTHREADS-SSE3 -s 70perc_raxml.phylip -n run2 -m GTRCATI -f a -N 100 -x $RANDOM -p $RANDOM -T 12
```
Following runs, downloaded RAxML_bipartitions file for each run and checked for topological convergence. After confirming this, run with the best likelihood (as presented in RAxML_info) used as representative tree.  

Installing exabayes
```
wget https://cme.h-its.org/exelixis/resource/download/software/exabayes-1.5.1.tar.gz
module load GCC/9.2.0
tar -zvxf exabayes-1.5.1.tar.gz
rm exabayes-1.5.1.tar.gz
```
Creating the config.nexus file exabayes needs (contents below)
```
begin run; 
   numRuns 4
   numCoupledChains 3
end;
```
Running exabayes (did these steps for both the 50 and 70% complete matrix - changing directories/input files etc etc)
```
#!/bin/bash -e

#SBATCH -A uoo00105 
#SBATCH -J 50perc_run1
#SBATCH --ntasks 1
#SBATCH -c 36
#SBATCH -t 48:00:00
#SBATCH --mem=105G
#SBATCH -D /nesi/nobackup/uoo00105/beetles/50perc_raxml
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread

module load GCC/9.2.0
/nesi/nobackup/uoo00105/beetles/exabayes-1.5.1/yggdrasil -f /nesi/nobackup/uoo00105/beetles/50perc_raxml/50perc_raxml.phylip -m DNA -s $RANDOM -n run1 -T 12 -M 0 -c config.nexus

```
If exabayes times out, you can restat using the -r flag in place of the -n flag
```
#!/bin/bash -e

#SBATCH -A uoo00105 
#SBATCH -J 50perc_run1
#SBATCH --ntasks 1
#SBATCH -c 36
#SBATCH -t 48:00:00
#SBATCH --mem=105G
#SBATCH -D /nesi/nobackup/uoo00105/beetles/50perc_raxml
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread

module load GCC/9.2.0
/nesi/nobackup/uoo00105/beetles/exabayes-1.5.1/yggdrasil -f /nesi/nobackup/uoo00105/beetles/50perc_raxml/50perc_raxml.phylip -m DNA -s $RANDOM -n run1_B -r run1 -T 12 -M 0 -c config.nexus
```
Summarising the parameter space for Exabayes using postProcParam
```
module load GCC/9.2.0
/nesi/nobackup/uoo00105/beetles/exabayes-1.5.1/postProcParam -n run1 -f ExaBayes_parameters*run1
/nesi/nobackup/uoo00105/beetles/exabayes-1.5.1/postProcParam -n run2 -f ExaBayes_parameters*run2

# Checking overlap of parameters between runs in R
library(tidyverse)

# Reading in file 
data <- read_tsv("../exabayes_convergence_results_5Jul2021.txt")

# Deleting superfluous read-in columns for 70% matrix (50% matrix fine - had deleted spaces before converting to tab delimited)
data <- data[,-c(13:17)]

data
## A tibble: 28 x 12
#   paramName       mean      sd     perc5    perc25    median    perc75     per95   eSS  psrf matrix   run
#   <chr>          <dbl>   <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl> <dbl> <dbl> <chr>  <dbl>
# 1 alpha{0}     2.44e-1 1.39e-3   2.44e-1   2.43e-1   2.44e-1   2.45e-1   2.46e-1  5380     1 70perc     1
# 2 r{0}(C<->…   3.51e-1 2.63e-3   3.51e-1   3.49e-1   3.51e-1   3.53e-1   3.55e-1  1710     1 70perc     1
# 3 r{0}(C<->…   8.28e-2 1.08e-3   8.28e-2   8.21e-2   8.28e-2   8.36e-2   8.46e-2  3910     1 70perc     1
# 4 r{0}(G<->…   6.71e-2 9.93e-4   6.71e-2   6.64e-2   6.71e-2   6.77e-2   6.87e-2  1510     1 70perc     1
# 5 TL{0}        6.75e+0 4.03e-2   6.75e+0   6.72e+0   6.75e+0   6.78e+0   6.82e+0  5290     1 70perc     1
# 6 pi{0}(G)     2.39e-1 1.41e-3   2.39e-1   2.38e-1   2.39e-1   2.4 e-1   2.42e-1  2490     1 70perc     1
# 7 r{0}(A<->…   7.59e-2 1.06e-3   7.59e-2   7.52e-2   7.59e-2   7.66e-2   7.76e-2  3430     1 70perc     1
# 8 pi{0}(A)     2.64e-1 1.46e-3   2.64e-1   2.63e-1   2.64e-1   2.65e-1   2.66e-1  2210     1 70perc     1
# 9 LnL         -6.89e+5 8.04e+0  -6.89e+5  -6.89e+5  -6.89e+5  -6.89e+5  -6.89e+5  3240     1 70perc     1
#10 r{0}(A<->…   1.03e-1 1.19e-3   1.03e-1   1.02e-1   1.03e-1   1.04e-1   1.05e-1  3520     1 70perc     1

for (i in unique(data$paramName)) {
  tempdata <- data[which(data$paramName==i),]
  ggplot(tempdata) +
    geom_boxplot(mapping=aes(x=as.factor(run), ymin=perc5, lower=perc25, middle=median, upper=perc75, max=per95),  stat = "identity") + ggtitle(i) + theme_bw() + theme(
      plot.title = element_text(hjust=0.5))
  ggsave(paste(i,".pdf"))
}
getwd()

```
Summarising the topologies for each run
```
module load GCC/9.2.0
/nesi/nobackup/uoo00105/beetles/exabayes-1.5.1/consense -n run1 -f ExaBayes_topologies*run1
/nesi/nobackup/uoo00105/beetles/exabayes-1.5.1/consense -n run2 -f ExaBayes_topologies*run2
```
Downloaded the newick tree files and checked for topological convergence, then merged tree files from each run:
```
module load GCC/9.2.0
/nesi/nobackup/uoo00105/beetles/exabayes-1.5.1/consense -n total -f ExaBayes_topologies*
```
