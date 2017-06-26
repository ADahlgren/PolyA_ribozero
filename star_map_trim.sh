module load star

ln -fs /share/finnolab/adahl/polyA_ribozero/qc_data/*_trim.qc.fq .

for filename in *R1_001_trim.qc.fq
do
base=$(basename $filename .qc.fq)
echo $base
base2=${base/R1/R2}
echo $base2

#map with star - outSAMstrandField includes tag XS
STAR --runThreadN 8 --genomeDir star_index \
--readFilesIn /share/finnolab/adahl/polyA_ribozero/mapping/${base}.qc.fq \
/share/finnolab/adahl/polyA_ribozero/mapping/${base2}.qc.fq \
--outFileNamePrefix alignments/${base}_2 --outSAMstrandField intronMotif --outSAMtype BAM Unsorted

done
