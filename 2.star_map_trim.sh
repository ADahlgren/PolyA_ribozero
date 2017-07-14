module load star

mkdir star_index

ln -fs /share/finnolab/adahl/polyA_ribozero/stringtie/*.GTF

#Genome indexing
#runMode genomeGenerate - telling star to "run genome indices generation job"
#genomeDir star_index - "specifies path to directory where genome indices are stored"
#genomeFastaFiles genome.fa - indicates the genome to use as reference
#sjdbGTFfile mergedTrans.GTF - indicates file with annotated transcripts in GTF
#sjdbOverhang 125 - According to manual, this should be ReadLength-1,in this case 126 - 1 = 125
STAR --runThreadN 8 --runMode genomeGenerate \
--genomeDir star_index --genomeFastaFiles genome.fa \
--sjdbGTFfile mergedTrans.GTF --sjdbOverhang 125


ln -fs /share/finnolab/adahl/polyA_ribozero/qc_data/*_trim.qc.fq .

mkdir alignments

for filename in *R1_001_trim.qc.fq
do
base=$(basename $filename .qc.fq)
echo $base
base2=${base/R1/R2}
echo $base2

#map with star
#outSAMstrandField includes tag XS needed for stringtie
#outSAMtype will be bam file sorted by coordinate since sorted file is needed for next steps
STAR --runThreadN 8 --genomeDir star_index \
--readFilesIn /share/finnolab/adahl/polyA_ribozero/mapping/${base}.qc.fq \
/share/finnolab/adahl/polyA_ribozero/mapping/${base2}.qc.fq \
--outFileNamePrefix alignments/${base}_ --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate

done
