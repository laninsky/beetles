The ```concatenation_script.R``` within this repo needs a folder of phylip files and a ```*uce_type_summary.txt``` generated using the python scripts from van Damn et al. (2018):  
Van Dam M. H., Henderson J. B., Esposito L., Trautwein M. 2018. Incorporating functional genomics into UCE phylogenomics improves gene and species tree reconstructions. Syst. Biol.  
(https://github.com/matthewhvandam/integrating-functional-genomics-into-phylogenomics). 

The ```concatenation_script.R``` uses the ```*uce_type_summary.txt``` file to concatenate loci that belong to the same gene and copies over the loci that are NOT concatenated to a subfolder, non_concat_loci.   

Assumptions: taxa are named the same way in all of the uce files, and this matches the way they are named in uce_summary_file_path. All taxa are present in all phylip files pointed to in phylip_w_missing_path (e.g. they are phylip files w missing data added). If there are files in uce_summary_file_path that aren't present in uce_summary_file_path the code should handle this OK.  

You'll need to modify the following lines of code at Line 9, 13, and 17 to point to the appropriate paths for your datasets:  

*Pathway to the *uce_type_summary.txt file*
```uce_summary_file_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/All_beetle_probeset_of_Faircloth/Coleoptera_Tricas_uce_type_summary.txt"```

*Pathway to the phylip input files*
```phylip_w_missing_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/All_beetle_probeset_of_Faircloth/50perc_internal_w_missing_phylip"```

*Where you would like the outfiles written to*
```outfile_path <- "/Users/alanaalexander/Dropbox/UCE_characterization_data_Alana/All_beetle_probeset_of_Faircloth/van_dam_concatenated"```

Some examples of paths are given here, and also in the code.

These scripts depend on the package tidyverse, and on base R.
