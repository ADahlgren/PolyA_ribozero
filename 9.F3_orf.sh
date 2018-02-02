#transdecoder to identify open reading frames (takes quite a bit of time)

module load transdecoder/3.0.0

for filename in *f1.fa
do
base=$(basename $filename .fa)
echo $base

#t = fasta file

TransDecoder.LongOrfs -t ${base}.fa

done

#output = directory ending in .transdecoder_dir; contains 5 files describing base freqs and longest orfs
