
module load stringtie

#link mapped bam files from star pass2 and gtf files from mapping file
ln -fs /share/finnolab/adahl/polyA_ribozero/mapping/pass2/*.out.bam
ln -fs /share/finnolab/adahl/polyA_ribozero/mapping/star_index/*.GTF

#create base
for filename in *.out.bam
do
base=$(basename $filename _Aligned.sortedByCoord.out.bam)
echo $base

#stringtie
echo "stringtie time"
stringtie ${base}_Aligned.sortedByCoord.out.bam -G mergedTrans.GTF -o ${base}_merged_stringtie.gtf -A ${base}_merged_abund.out \
-b /share/finnolab/adahl/polyA_ribozero/ballgown/merged_gtf/${base}_ballgown

done
