for filename in *tab.bed
do
base=$(basename $filename _tab.bed)
echo $base

#filter out all transcripts with TPM < 0.1 and with length < 199
#filter based on length - >=200 based on definition of lncRNA (note: there were quite a few transcripts at 199bp)

cat ${base}_tab.bed | awk '{if ($6 >= 0.1) {print $0}}' | awk '{if ($7 >= 200) {print $0}}' > ${base}_filter12.bed

done
