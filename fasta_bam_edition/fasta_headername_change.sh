# ths script replaces the header of the fasta sequences  
#name 


for file in *.fas;
do 
sed -i "s/^>.*/>"${file%.*}"/gi" "$file";
done
