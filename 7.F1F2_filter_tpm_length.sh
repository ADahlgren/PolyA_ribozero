for filename in *tab.bed
do
base=$(basename $filename _tab.bed)
echo $base

#filter out all transcripts with TPM < 0.1 and with length < 199
#filter based on length - >=199 since lncRNA are longer than that size and a bunch of transcripts were 199 bp

cat ${base}_tab.bed | awk '{if ($6 >= 0.1) {print $0}}' | awk '{if ($7 >= 199) {print $0}}' > ${base}_filter12.bed

done
