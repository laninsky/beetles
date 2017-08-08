fastqlineno=`wc -l $filename | awk '{print $1}'`
if [ $fastqlineno > 2147483647 ];
then echo "You will need to split your fastq file in two to proceed as it has too many lines to be read in by R"; 
else sampleevery=`echo "scale=0 ; 1 / $perc_to_sample" | bc`;
seq 1 $(( sampleevery * 4 )) $fastqlineno > linestosample;
echo $filename > filename;
Rscript downsampling_fastq.R;
fi
