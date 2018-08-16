## Following scripts from E. Scott and T. Mansour https://github.com/eyscott/lncRNA


for filename in *_sick_Aligned.sortedByCoord_merged_stringtie.gtf #output from stringtie
do
base=$(basename $filename _R1_001_sick_Aligned.sortedByCoord_merged_stringtie.gtf)
echo $base

#separate transcript and exon rows
awk '$3 ~ /transcript/' ${base}_R1_001_sick_Aligned.sortedByCoord_merged_stringtie.gtf > ${base}_transcripts.gtf

#extract column with a lot of info including transcript ids
#isolate two columns that may contain transcript id and manipulate format to allow for extraction of info
#if column 3=1, column 4 has transcript id, if not column 2 has transcript id
awk 'BEGIN {FS="\t"; OFS="\t";} {print $9}' ${base}_transcripts.gtf | \
awk 'BEGIN {FS=";"; OFS=" ";} {print $2" "$3}' | sed 's/reference_id/1/g' | \
sed 's/"//g' | sed 's/cov/2/g' | sed 's/ /\t/g' | \
awk '{if($3==1){print $4}else{print $2}}' > ${base}_transcript_ids.txt

#extract column with a lot of info including tpms
#isolate two columns that may have tmp and manipulate format to allow for easy extraction of info
#if column 1=2, tmp is in column 4, if not, tpm is in column 2
awk 'BEGIN {FS="\t"; OFS="\t";} {print $9}' ${base}_transcripts.gtf | \
awk 'BEGIN {FS=";"; OFS=" ";} {print $5" "$7}' | \
sed 's/TPM/1/g' | sed 's/cov/2/g' | sed 's/"//g' | sed 's/ /\t/g' | \
awk '{if($1==2){print $4}else{print $2}}' > ${base}_tpms.txt

awk '{print $7}' ${base}_transcripts.gtf > ${base}_strand.bed

awk '{print $1"\t"$4"\t"$5}' ${base}_transcripts.gtf > ${base}_bed1.bed

#paste together all necessary info in specific order
paste ${base}_bed1.bed ${base}_transcript_ids.txt ${base}_strand.bed ${base}_tpms.txt > ${base}.bed

wait

#calculate length
awk '{$7 = $3 - $2} 1' ${base}.bed > ${base}_length.bed

#remove header and make tab deliminated
sed '1d;s/ /\t/g' ${base}_length.bed > ${base}_tab.bed

done
