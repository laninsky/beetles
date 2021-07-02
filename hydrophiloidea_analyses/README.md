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
