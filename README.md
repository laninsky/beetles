# downsampling_fastq

arguments to export:
```
export filename=/panfs/pfs.local/scratch/bi/trunks18/For_Cluster/reads/full_bemHap/bemHap1-pe100-reads.fq

export perc_to_sample=0.25
```

If you get the following message:
"You will need to split your fastq file in two to proceed as it has too many lines to be read in by R"

You'll need to use split_downsampling_fastq.sh instead of downsampling_fastq.sh. The Rscript should remain the same.
