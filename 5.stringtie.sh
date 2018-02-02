#note: this is for refined GTF which I used in lnc analysis

module load stringtie

#link mapped bam files from star pass2 and gtf files from mapping file
ln -fs /share/finnolab/adahl/polyA_ribozero/mapping/pass2/*.out.bam
ln -fs /share/finnolab/adahl/polyA_ribozero/mapping/star_index/*.GTF

for filename in *.bam
do
base=$(basename $filename .out.bam)
echo $base

#already make ballgown directory and then this is to make subdirectories separating ballgown files for each sample
#may use these files later to look at differential expression?
cd ballgowns
mkdir ${base}_ballgown
cd ../

#stringtie
#G - GTF file (used both merged and filtered, proceded with filtered to lnc analysis)
#o - output gtf name
#A - report gene abundances in abund.out file
#b - path to where ballgown files should be outputed 
echo "stringtie time"
stringtie ${base}.out.bam -G filtered4_refined_Alltissues.GTF -o ${base}_refined_stringtie.gtf -A ${base}_refined_abund.out \
-b /share/finnolab/adahl/polyA_ribozero/ballgown/refined2/ballgowns/${base}_ballgown

done
