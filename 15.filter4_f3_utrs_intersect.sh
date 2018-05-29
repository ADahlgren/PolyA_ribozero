#There were quotes around some of the columns, and this was the easiest way to get rid of them

for filename in *_5.bed
do
base=$(basename $filename .bed)
echo $base

sed 's/\"//g' ${base}.bed > ${base}_unquoted.bed

done


for filename in *_3.bed
do
base=$(basename $filename .bed)
echo $base

sed 's/\"//g' ${base}.bed > ${base}_unquoted.bed

done

############################################################

#use bedtools intersect to identify overlap between output of filter3 and UTRs
# -c = report number of overlaps file a has with file b 
# -a = file a, compared to file b
# -b = designates file b
# -sorted = requires position-sorted input, helps keep memory usage low
# -s = force strandedness, maybe I shouldn't do this as there are many lines that have a . instead of +/-

module load bedtools2/2.27.0

for filename in *f3.bed
do
base=$(basename $filename _f3.bed)
echo $base

bedtools intersect -c -a ${base}_f3.bed -b ${base}_refined_nolncRNA_5.bed -sorted -s > ${base}_5.bed

bedtools intersect -c -a ${base}_f3.bed -b ${base}_refined_nolncRNA_3.bed -sorted -s > ${base}_3.bed

done
