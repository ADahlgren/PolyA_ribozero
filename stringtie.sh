module load samtools
module load stringtie

for filename in *_2Aligned.out.bam
do
base=$(basename $filename _2Aligned.out.bam)
echo $base

#sort bam files
echo "sorting"
samtools sort ${base}_2Aligned.out.bam > ${base}_sorted.bam

wait

#stringtie
echo "stringtie time"
stringtie ${base}_sorted.bam -G filtered4_refined_Alltissues.GTF > ${base}_refined_stringtie.gtf
stringtie ${base}_sorted.bam -G mergedTrans.GTF > ${base}_merged_stringtie.gtf

done
