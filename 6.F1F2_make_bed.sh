## Following scripts from E. Scott and T. Mansour https://github.com/eyscott/lncRNA


for filename in *sick*.out
do
base=$(basename $filename _R1_001_sick_Aligned.sortedByCoord_merged_abund.out)
echo $base

#extract only columns I want
#chr, start position, end position, gene ID, strand, TPM
cat ${base}_R1_001_sick_Aligned.sortedByCoord_merged_abund.out | awk 'BEGIN {FS="\t"; OFS="\t";} \
{print $3, $5, $6, $1, $4, $9}' > ${base}.bed

#calculate length
awk '{$7 = $3 - $2} 1' ${base}.bed > ${base}_length.bed

#remove header and make tab deliminated
sed '1d;s/ /\t/g' ${base}_length.bed > ${base}_tab.bed

done
