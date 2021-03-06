module load star/2.5.2b

mkdir star_index

#link reference files to index folder
cd star_index

#GTF file was downloaded from horse_trans github

#Genome indexing
  #runMode genomeGenerate - telling star to "run genome indices generation job"
  #genomeDir star_index - "specifies path to directory where genome indices are stored"
  #genomeFastaFiles genome.fa - indicates the genome to use as reference
  #sjdbGTFfile mergedTrans.GTF - indicates file with annotated transcripts in GTF
  #sjdbOverhang 125 - According to manual, this should be ReadLength-1,in this case 126 - 1 = 125
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star_index \
--genomeFastaFiles star_index/genome_3.0.fna --sjdbGTFfile star_index/GCF_002863925.1_EquCab3.0_genomic.gff \
--sjdbOverhang 125
