fastqlineno=`wc -l $filename | awk '{print $1}'`

split -l $(( fastqlineno / 2)) $filename

for i in xa*;
do fastqlineno=`wc -l $i | awk '{print $1}'`
if [ $fastqlineno > 2147483647 ];
then echo "You will need to split your fastq file further in order to proceed as it still has too many lines to be read in by R"; 
else sampleevery=`echo "scale=0 ; 1 / $perc_to_sample" | bc`;
seq 1 $(( sampleevery * 4 )) $fastqlineno > linestosample;
echo $i > filename;
Rscript downsampling_fastq.R;
fi;
done;
