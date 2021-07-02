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
conda activate /nesi/nobackup/uoo00105/beetles/conda
phyluce_align_get_only_loci_with_min_taxa --alignments abyss_50perc_nexus --taxa 63 --percent 0.7 --output 70perc_nexus --cores 12 --log-path logs
```
