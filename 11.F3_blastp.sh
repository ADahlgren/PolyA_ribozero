#run blastp
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
#UniProt is central access point for curated protein infor, including function, classification, and cross-references
#knowlegebase complete consists of fully annoted entries and computer-generated entries with automated classifcation and annotation in fasta format
module load blast/2.6.0+
#produces BLAST databases from fasta files

makeblastdb -in uniprot_sprot.fasta -dbtype prot

#sequence identification and similarity searches
#query = input file name
#dp = BLAST data base name
#max_target_seqs = maxiumum number of aligned sequences to keep
#outfmt = format of output (see http://www.metagenomics.wiki/tools/blast/blastn-output-format-6)
#evalue = expectation value threshold for saving hits
blastp -query 683610_Liver_Ribozero_S51_f1.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > 683610_Liver_Ribozero_sprot.outfmt6
blastp -query 683610_Liver_S61_f1.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > 683610_Liver_sprot.outfmt6
blastp -query 683610_Parietal_Cortex_Ribozero_S53_f1.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > 683610_Parietal_Cortex_Ribozero_sprot.outfmt6
blastp -query 683610_Parietal_Cortex_S55_f1.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > 683610_Parietal_Cortex_sprot.outfmt6
blastp -query 686521_Liver_Ribozero_S50_f1.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > 686521_Liver_Ribozero_sprot.outfmt6
blastp -query 686521_Liver_S60_f1.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > 686521_Liver_sprot.outfmt6
blastp -query 686521_Parietal_Cortex_Ribozero_S52_f1.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > 686521_Parietal_Cortex_Ribozero_sprot.outfmt6
blastp -query 686521_Parietal_Cortex_S54_f1.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 >  686521_Parietal_Cortex_sprot.outfmt6
