#download Pfam-A db from:  ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/
#Pfam is a database of protein domain models
curl -O ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm.gz
#unzip file
gunzip Pfam-A.hmm.gz

#Predict likely coding regions using ORFs predicted by transdecoder
#hmmsearch: search protein profile against protein sequence database


module load hmmer/3.1b2

hmmsearch --cpu 8 --tblout 683610_Liver_Ribozero_pfam.tblout Pfam-A.hmm 683610_Liver_Ribozero_S51_f1.fa.transdecoder_dir/longest_orfs.pep

hmmsearch --cpu 8 --tblout 683610_Liver_pfam.tblout Pfam-A.hmm 683610_Liver_S61_f1.fa.transdecoder_dir/longest_orfs.pep

hmmsearch --cpu 8 --tblout 683610_Parietal_Cortex_Ribozero_pfam.tblout Pfam-A.hmm 683610_Parietal_Cortex_Ribozero_S53_f1.fa.transdecoder_dir/longest_orfs.pep

hmmsearch --cpu 8 --tblout 683610_Parietal_Cortex_pfam.tblout Pfam-A.hmm 683610_Parietal_Cortex_S55_f1.fa.transdecoder_dir/longest_orfs.pep

hmmsearch --cpu 8 --tblout 686521_Liver_Ribozero_pfam.tblout Pfam-A.hmm 686521_Liver_Ribozero_S50_f1.fa.transdecoder_dir/longest_orfs.pep

hmmsearch --cpu 8 --tblout 686521_Liver_pfam.tblout Pfam-A.hmm 686521_Liver_S60_f1.fa.transdecoder_dir/longest_orfs.pep

hmmsearch --cpu 8 --tblout 686521_Parietal_Cortex_Ribozero_pfam.tblout Pfam-A.hmm 686521_Parietal_Cortex_Ribozero_S52_f1.fa.transdecoder_dir/longest_orfs.pep

hmmsearch --cpu 8 --tblout 686521_Parietal_Cortex_pfam.tblout Pfam-A.hmm 686521_Parietal_Cortex_S54_f1.fa.transdecoder_dir/longest_orfs.pep
