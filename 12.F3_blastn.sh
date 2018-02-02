#compare to human transcriptome; remove all that appear as cDNA
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

#isolate positions that are protein coding in human 
grep "gene_biotype:protein_coding" Homo_sapiens.GRCh38.cdna.all.fa | awk -F"[> ]" '{print $2;}' > protein_coding.ids

module load qiime/1.9.1
#filter_fasta = remove sequences from a fasta or fastq file based on input criteria
#f - input fasta
#o - output fasta
#s - A list of sequence identifiers (or tab-delimited lines with a seq identifier in the first field) which should be retained
filter_fasta.py -f Homo_sapiens.GRCh38.cdna.all.fa -o Homo_sapiens.GRCh38.cdna.ptn.fa -s protein_coding.ids

################################################################################################################3

module load blast/2.6.0+

#create BLAST databases
#in - input file
#dbtype - specifiy database (nucleotide)
#out - output database name
echo "make blast db -> hg38_cdnaPtn_db"
makeblastdb -in Homo_sapiens.GRCh38.cdna.ptn.fa -dbtype nucl -out hg38_cdnaPtn_db

wait

for filename in *f1.fa
do
base=$(basename $filename .fa)
echo ${base}

#query - input file
#db - database name
#max_target_seqs - kind of self explanitory, followed by number of sequences
#outfmt - output format (see http://www.metagenomics.wiki/tools/blast/blastn-output-format-6)
#evalue - see http://www.metagenomics.wiki/tools/blast/evalue
#num_treads - number of threads obviously

blastn -query ${base}.fa  -db hg38_cdnaPtn_db  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 4 > ${base}_hg38_cdna.outfmt6

done
