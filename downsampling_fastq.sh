fastqlineno=`wc -l $filename | awk '{print $1}'`
echo "Your file has the following number of lines" $fastqlineno
if (( $fastqlineno > 2147483647 ));
then
echo "You will need to split your fastq file in two to proceed";
echo "as it has too many lines to be read in by R. It has the following no of lines";
echo $fastqlineno;
else
sampleevery=`echo "scale=0 ; 1 / $perc_to_sample" | bc`;
echo "We will sample every " $sampleevery " sequences";
seq 1 $(( sampleevery * 4 )) $fastqlineno > linestosample;
echo $filename > filename;
Rscript downsampling_fastq.R;
fi;
