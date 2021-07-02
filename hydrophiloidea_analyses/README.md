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
conda init bash
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
Running exabayes
```
#!/bin/bash -e

#SBATCH -A uoo00105 
#SBATCH -J 50perc_run1
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 12:00:00
#SBATCH --mem=20G
#SBATCH -D /nesi/nobackup/uoo00105/beetles/50perc_raxml
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread

/nesi/nobackup/uoo00105/beetles/exabayes-1.5.1/yggdrasil -f /nesi/nobackup/uoo00105/beetles/50perc_raxml/50perc_raxml.phylip -m DNA -s $RANDOM -n run1 -T 12 -M 0 -c config.nexus

```
If exabayes times out, you can restat using the -r flag in place of the -n flag
```
#!/bin/bash -e

#SBATCH -A uoo00105 
#SBATCH -J 50perc_run1
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 12:00:00
#SBATCH --mem=20G
#SBATCH -D /nesi/nobackup/uoo00105/beetles/50perc_raxml
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread

/nesi/nobackup/uoo00105/beetles/exabayes-1.5.1/yggdrasil -f /nesi/nobackup/uoo00105/beetles/50perc_raxml/50perc_raxml.phylip -m DNA -s $RANDOM -n run1_B -r run1 -T 12 -M 0 -c config.nexus
```
