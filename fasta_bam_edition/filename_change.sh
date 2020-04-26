# This scripts modifies the name of the files if they have numbers and underscore before the ID of the samples.
# Also it changes .20 and .2, which corresponds to 20 and 2 percentage  to _20PC and _2PC

shopt -s extglob


for fname in  *.fas;
do mv -- "$fname" "${fname##*([[:digit:]])_}";
done 

for file in  *.fas;
do 

mv "$file" "${file//.20.fas/_20PC.fas}";
mv "$file" "${file//.2.fas/_2PC.fas}";

done
