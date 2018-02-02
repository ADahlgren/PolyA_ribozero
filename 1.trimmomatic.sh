#load trimmomatic
module load trimmomatic/0.33
#load fastqc
module load fastqc/v0.11.2

for filename in *R1_001.fastq.gz
do
#make base by removing R1_001.fastq.gz
base=$(basename $filename .fastq.gz)
echo $base
#construct R2 filename
base2=${base/R1/R2}
echo $base2

#run trimmomatic - PE for paired end; 
trimmomatic PE ${base}.fastq.gz ${base2}.fastq.gz ${base}_trim.qc.fq ${base}_trim_s1_se ${base2}_trim.qc.fq ${base2}_trim_s2_se SLIDINGWINDOW:3:28

#fastqc
fastqc ${base}_trim.qc.fq
fastqc ${base2}_trim.qc.fq

done
