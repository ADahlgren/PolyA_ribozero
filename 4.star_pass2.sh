#COPY TAB FILES (from pass1) INTO THIS

module load star

mkdir pass2

for filename in 683610_Liver*R1_001_trim.qc.fq
do
base=$(basename $filename .qc.fq)
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
STAR --runThreadN 12 --genomeDir star_index \
--sjdbFileChrStartEnd /share/finnolab/adahl/polyA_ribozero/mapping/pass1/683610_Liver_Ribozero_S51_R1_001_trim_SJ.out.tab /share/finnolab/adahl/polyA_ribozero/mapping/pass1/683610_Liver_S61_R1_001_trim_SJ.out.tab \
--readFilesIn /share/finnolab/adahl/polyA_ribozero/mapping/${base}.qc.fq \
/share/finnolab/adahl/polyA_ribozero/mapping/${base2}.qc.fq \
--outFileNamePrefix pass2/${base}_ --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate

done
