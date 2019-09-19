#note: this is for refined GTF which I used in lnc analysis

module load stringtie/1.3.2d

#link mapped bam files from star pass2 and gtf files from mapping file

for filename in *.bam
do
base=$(basename $filename .out.bam)
echo $base


#stringtie
#G - GTF file (used both merged and filtered, proceded with filtered to lnc analysis)
#o - output gtf name
#A - report gene abundances in abund.out file
 
echo "stringtie time"
stringtie ${base}.out.bam -G ref_EquCab3.0_top_level.gtf.GeneBank -o ${base}_equcab3_stringtie.gtf -A ${base}_equcab3_abund.out

done
