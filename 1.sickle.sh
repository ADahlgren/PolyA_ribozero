#load sickle
module load sickle
#load fastqc
module load fastqc

for filename in *R1_001.fastq.gz
do
#make base by removing R1_001.fastq.gz
base=$(basename $filename .fastq.gz)
echo $base
#construct R2 filename
base2=${base/R1/R2}
echo $base2

#run sickle
#f - forward read
#r - reverse read
#t - sanger even though illumina was actually used

sickle pe -f ${base}.fastq.gz -r ${base2}.fastq.gz -t sanger -o ${base}_sick.fastq -p ${base2}_sick.fastq -s ${base}_sick_singles.fastq

#fastqc
fastqc ${base}_sick.fastq
fastqc ${base2}_sick.fastq

done
