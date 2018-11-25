for filename in *length.bed
do
base=$(basename $filename _length.bed)
echo $base

#filter out all transcripts with TPM < 0.5 (based on histogram) and with length < 200 (based on definition of lncRNA
#though there are quite a few at 199 bp
cat ${base}_length.bed | awk '{if ($6 >= 0.5) {print $0}}' | awk '{if ($7 >= 200) {print $0}}' > ${base}_filter12.bed

done
