# generate_downsampled v0.0
The Rscript in this repository, generate_downsampled.R, contains a function that can be used to output downsampled datasets of fasta files (i.e. bootstrapping without replacement). To use, either download and source the script in your R console, or copy the function into your R console.

As arguments it expects:
* the path to a folder containing fasta files, 'fasta_location'
* The location where you want to put the downsampled datasets, 'output_dir'
* The number of replicate downsampled datasets you want to create, 'reps'
* The number of loci you want to downsample to, 'loci_per_rep'  
e.g.  
generate_downsampled(fasta_location,output_dir,reps,loci_per_rep)

Example usage:
```
generate_downsampled("/Users/alanaalexander/Dropbox/UCE_data_for_Alana/50perc_data_matrix/50perc_internal_fasta","/Users/alanaalexander/Dropbox/beetles/75_downsamples",100,305)
```
This would create 100 replicate folders within `/Users/alanaalexander/Dropbox/beetles/75_downsamples` each containing 305 fasta files sampled from the total number of fasta files available at `/Users/alanaalexander/Dropbox/UCE_data_for_Alana/50perc_data_matrix/50perc_internal_fasta`

# Version history  
v0.0: the code used in the following, where generate_downsampled was first published:
Gustafson, G.T., Baca, S.M., Alexander, A.M. and Short, A.E., 2020. Phylogenomic analysis of the beetle suborder Adephaga with comparison of tailored and generalized ultraconserved element probe performance. Systematic Entomology, 45(3): 552-570.

If you want to cite the specific version of the script you used, I suggest the following as well as the pub above:  
Alexander, A. 2019. downsampled_dataset vx.x. Available from https://github.com/laninsky/beetles/new/master/downsampled_dataset

This script wouldn't be possible without:  
R: R Core Team. 2015. R: A language and environment for statistical computing. URL http://www.R-project.org/. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/

