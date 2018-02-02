module load bedtools2/2.27.0

for filename in *filter12.bed
do
base=$(basename $filename _filter12.bed)
echo $base


#make fasta file from bed (you get sequence for each region in bed file)
#s = force strandedness, so if feature is on antisense strand, it'll be reverse complemented
#fi = precedes reference fasta
#bed = precedes bed file
#name = "Use the “name” column in the BED file for the FASTA headers in the output FASTA file"
#fo = output file name

bedtools getfasta -s -fi genome.fa -bed ${base}_filter12.bed \
-name -fo ${base}_f1.fa

done
