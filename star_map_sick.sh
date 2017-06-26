module load star

ln -fs /share/finnolab/adahl/polyA_ribozero/qc_data/*_sick.fastq .

for filename in *R1_001_sick.fastq
do
base=$(basename $filename .fastq)
echo $base
base2=${base/R1/R2}
echo $base2

#map with STAR
STAR --runThreadN 16 --genomeDir star_index \
--readFilesIn /share/finnolab/adahl/polyA_ribozero/mapping/${base}.fastq \
/share/finnolab/adahl/polyA_ribozero/mapping/${base2}.fastq \
--outFileNamePrefix alignments/${base}_2 --outSAMstrandField intronMotif --outSAMtype BAM Unsorted

done
