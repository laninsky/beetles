fastqlineno=`wc -l $filename | awk '{print $1}'`
sampleevery=`echo "scale=0 ; 1 / $perc_to_sample" | bc`
seq 1 $(( sampleevery * 4 )) $fastqlineno > linestosample
echo $filename > filename
Rscript downsampling_fastq.R
