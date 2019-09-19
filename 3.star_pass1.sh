module load star/2.5.2b

mkdir pass1

#make base (this is just subset for now)
for filename in *R1_001_*.fastq
do
base=$(basename $filename .fastq)
echo $base
base2=${base/R1/R2}
echo $base2

#Map with star (first pass)
  #genomeDir - location of reference
  #readFilesIn - path to each read
  #outFileNamePrefix - location and name of output
  #outSAMstrandField intronMotif - for XS strand attribute needed for stringtie
  #outSAMtype - type of output; unsorted bam file (will sort during second pass)

STAR --runThreadN 10 --genomeDir /share/finnolab/reference/equCab3.0/NCBI/indcies/star_index \
--readFilesIn ${base}.fastq ${base2}.fastq \
--outFileNamePrefix pass1/${base}_ --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate
done

#copy .tab files into pass2
