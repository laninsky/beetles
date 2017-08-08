fastqlineno=`wc -l $filename | awk '{print $1}'`
sampleevery=`$(( 1 / perc_to_sample ))`
sampleevery=`$(( sampleevery * 4 ))`
seq 1 $sampleevery $fastqlineno > linestosample
echo $filename > filename
Rscript downsampling_fastq.R
