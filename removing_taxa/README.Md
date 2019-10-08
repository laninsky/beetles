# removing_phylip_taxon v0.0  
The Rscript in this repository, removing_phylip_taxon.R, contains a function that can be used to output phylip files after removing a specific taxon. To use, either download and source the script in your R console, or copy the function into your R console.  

As arguments it expects:  
the path to a folder containing phylip files, 'working_dir'  
the taxon you want to remove, 'taxon_name'  
e.g.  
`removing_phylip_taxon(working_dir,taxon_name) ` 

Example usage:  
`removing_phylip_taxon("/Users/alanaalexander/Dropbox/UCE_data_for_Alana/Spanglerogyrus_removal/50perc_phylip_alignments_w_Spangy","sle1702_sp")`
This would create a folder "sle1702_sp_removed" in the directory the Rscript was launched from, containing phylip files as at /Users/alanaalexander/Dropbox/UCE_data_for_Alana/Spanglerogyrus_removal/50perc_phylip_alignments_w_Spangy, except non-interleaved, and with the taxon sle1702_sp removed (and the number of taxa in the phylip heaer line adjusted).  

## Version history
v0.0: the code used in Beetles TBD, where generate_downsampled was first published.

If you want to cite the specific version of the script you used, I suggest the following as well as the pub above:
Alexander, A. 2019. removing_phylip_taxon vx.x. Available from https://github.com/laninsky/beetles/new/master/removing_taxa

This script wouldn't be possible without:
R: R Core Team. 2015. R: A language and environment for statistical computing. URL http://www.R-project.org/. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/
