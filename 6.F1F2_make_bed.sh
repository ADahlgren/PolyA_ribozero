## Following scripts from E. Scott and T. Mansour https://github.com/eyscott/lncRNA


for filename in *_R1_001_sick_Aligned.sortedByCoord_equcab3_stringtie.gtf
do
base=$(basename $filename _R1_001_sick_Aligned.sortedByCoord_equcab3_stringtie.gtf)
echo $base

awk '$3 ~ /transcript/' ${base}_R1_001_sick_Aligned.sortedByCoord_equcab3_stringtie.gtf > ${base}_transcripts.gtf

awk 'BEGIN {FS="\t"; OFS="\t";} {print $9}' ${base}_transcripts.gtf | \
awk 'BEGIN {FS=";"; OFS=" ";} {print $2" "$3}' | sed 's/reference_id/1/g' | \
sed 's/"//g' | sed 's/cov/2/g' | sed 's/ /\t/g' | \
awk '{if($3==1){print $4}else{print $2}}' > ${base}_transcript_ids.txt

awk 'BEGIN {FS="\t"; OFS="\t";} {print $9}' ${base}_transcripts.gtf | \
awk 'BEGIN {FS=";"; OFS=" ";} {print $5" "$8}' | sed 's/TPM/1/g' | \
sed 's/ref_gene_name/2/g' |  sed 's/"//g' | sed 's/ /\t/g' | \
awk '{if($3==1){print $4}else{print $2}}' > ${base}_tpms.txt

awk '{print $7}' ${base}_transcripts.gtf > ${base}_strand.bed

awk '{print $1"\t"$4"\t"$5}' ${base}_transcripts.gtf > ${base}_bed1.bed


paste ${base}_bed1.bed ${base}_transcript_ids.txt ${base}_strand.bed ${base}_tpms.txt > ${base}.bed

done

