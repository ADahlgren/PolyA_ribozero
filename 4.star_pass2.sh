#COPY TAB FILES (from pass1) INTO THIS

module load star/2.5.2b

mkdir pass2

for filename in *R1_001_sick.fastq
do
base=$(basename $filename .fastq)
echo $base
base2=${base/R1/R2}
echo $base2

#COPY TAB FILES!

#map with star(pass2)
  #genomeDir - location of reference genome
  #sjdbFileChrStartEnd - tab files from pass1 (use junctions detected in 1st pass as ”annotated” junctions for 2nd pass)
  #readFilesIn - path to reads
  #outFileNamePrefix - location and name of output
  #outSAMstrandField intronMotif - for XS strand attribute needed for stringtie
  #outSAMtype - type of output; sorted bam file (so don't have to samtools sort later)
STAR --runThreadN 10 --genomeDir /share/finnolab/reference/equCab3.0/NCBI/indcies/star_index \
--sjdbFileChrStartEnd pass1/683610_Liver_Ribozero_S51_R1_001_sick_SJ.out.tab \
pass1/683610_Parietal_Cortex_S55_R1_001_sick_SJ.out.tab \
pass1/686521_Parietal_Cortex_Ribozero_S52_R1_001_sick_SJ.out.tab \
pass1/683610_Liver_S61_R1_001_sick_SJ.out.tab \
pass1/686521_Liver_Ribozero_S50_R1_001_sick_SJ.out.tab \
pass1/686521_Parietal_Cortex_S54_R1_001_sick_SJ.out.tab \
pass1/683610_Parietal_Cortex_Ribozero_S53_R1_001_sick_SJ.out.tab \
pass1/686521_Liver_S60_R1_001_sick_SJ.out.tab \
--readFilesIn ${base}.fastq ${base2}.fastq \
--outFileNamePrefix pass2/${base}_ --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate
done
