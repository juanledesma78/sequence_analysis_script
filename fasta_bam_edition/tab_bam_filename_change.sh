# This scripts modifies the name of the bam and tabular files 
#if they have numbers and underscore before the ID of the samples.

shopt -s extglob


for fname in  *.bam;
do mv -- "$fname" "${fname##*([[:digit:]])_}";

done

for fname in  *.tabular;
do mv -- "$fname" "${fname##*([[:digit:]])_}";

done
