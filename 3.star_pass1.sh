module load star/2.5.2b

mkdir pass1

#make base (this is just subset for now)
for filename in 683610_Liver*R1_001_trim.qc.fq
do
base=$(basename $filename .qc.fq)
echo $base
base2=${base/R1/R2}
echo $base2

#Map with star (first pass)
  #genomeDir - location of reference
  #readFilesIn - path to each read
  #outFileNamePrefix - location and name of output
  #outSAMstrandField intronMotif - for XS strand attribute needed for stringtie
  #outSAMtype - type of output; unsorted bam file (will sort during second pass)

STAR --runThreadN 12 --genomeDir star_index \
--readFilesIn /share/finnolab/adahl/polyA_ribozero/mapping/${base}.qc.fq \
/share/finnolab/adahl/polyA_ribozero/mapping/${base2}.qc.fq \
--outFileNamePrefix pass1/${base}_ --outSAMstrandField intronMotif --outSAMtype BAM Unsorted
done

#copy .tab files into pass2
