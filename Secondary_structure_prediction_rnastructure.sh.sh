#!/bin/bash

#SBATCH --job-name=Fold_Species
#SBATCH --partition=cpu6348
#SBATCH -n 5
#SBATCH --qos cpu6348
#SBATCH --mail-type=ALL
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --exclusive


ml load rnastructure
ml load viennarna
ml load parallel
ml load r

export DATAPATH=/gpfs/spack/opt/linux-rocky8-icelake/gcc-8.5.0/rnastructure-6.4-yyu25welaovp5esuz6jafgl5fsi3zvct/data_tables
#######################The first part is for unmodified sequence 
input_file=""mat_Species_sequence.fa"" ###This is an example input
output_dir="RNAstructure_fasta"

mkdir -p $output_dir 

awk -v RS=">" -v ORS="" '{if (NF) {print ">"$0 > "'$output_dir'/sequence_"++count".fasta"}}' $input_file

echo "FASTA file split into individual sequence files in the '$output_dir' directory."

cd RNAstructure_fasta

ls *.fasta > fastaname1.txt
cat fastaname1.txt | awk -F"." '{print $1}' > fastaname.txt
cat fastaname.txt | parallel -j 5 Fold {}.fasta {}.ct -mfe 
cat fastaname.txt | parallel -j 5 ct2dot {}.ct 1 {}.fold
cat fastaname.txt| parallel -j 5 cat {}.fold >> combined_seq
cat fastaname.txt | parallel -j 5 rm {}.fold
cat fastaname.txt | parallel -j 5 rm {}.fasta

mkdir ../CTa
cat fastaname.txt| parallel -j 5 mv {}.ct ../CTa
rm fastaname1.txt
rm fastaname.txt

#rename
input_file="combined_seq"
awk 'NR % 3 == 1 {print substr($0, 2)}' $input_file > MFE_structure.txt
awk 'NR % 3 == 1 {sub(".*m6A", "m6A", $0); print} NR % 3 != 1 {print}' $input_file > RNAstructure
awk 'NR % 3 == 1 {$0 = ">"$0"S"} {print}' RNAstructure > RNAstructureS

rm RNAstructure
rm combined_seq

mkdir ../combine
mv RNAstructureS ../combine
cd ..

#########################################################
#######################The second part is for modified sequence 
input_file="mat_Species_sequence_modification.fa" ###This is an example input
output_dir="RNAstructure_fastam"

mkdir -p $output_dir

awk -v RS=">" -v ORS="" '{if (NF) {print ">"$0 > "'$output_dir'/sequence_"++count".fasta"}}' $input_file

echo "FASTA file split into individual sequence files in the '$output_dir' directory."

cd RNAstructure_fastam
ls *.fasta > fastaname1.txt
cat fastaname1.txt | awk -F"." '{print $1}' > fastaname.txt
cat fastaname.txt | parallel -j 5 Fold {}.fasta {}.ct -mfe --alphabet m6A
cat fastaname.txt | parallel -j 5 ct2dot {}.ct 1 {}.fold
cat fastaname.txt| parallel -j 5 cat {}.fold >> combined_seq1
cat fastaname.txt | parallel -j 5 rm {}.fold
cat fastaname.txt | parallel -j 5 rm {}.fasta

mkdir ../CTm
cat fastaname.txt| parallel -j 5 mv {}.ct ../CTm
rm fastaname1.txt
rm fastaname.txt

#rename
input_file="combined_seq1"

awk 'NR % 3 == 1 {print substr($0, 2)}' $input_file > MFE_structure_m6A.txt
awk 'NR % 3 == 1 {sub(".*m6A", "m6A", $0); print} NR % 3 != 1 {print}' $input_file > RNAstructure_m6A1
awk 'NR % 3 == 1 {$0 = ">"$0"M"} {print}' RNAstructure_m6A1 > RNAstructure_m6A2

rm RNAstructure_m6A1
awk 'NR%3==2 {gsub("M", "A")} 1' RNAstructure_m6A2 > RNAstructure_m6A
rm RNAstructure_m6A2
rm combined_seq1

mv RNAstructure_m6A ../combine

cd ../combine
cat RNAstructureS RNAstructure_m6A > combined.fold

rm RNAstructureS
rm RNAstructure_m6A


awk '/^>/{if(x){close(x)} x=$0; sub("^>", "", x);}{print > x;}' combined.fold
rm *.fold
cat * > the_file_for_foster.fold
mv the_file_for_foster.fold ..
rm *

cd ..
rmdir combine

ml load viennarna
ml load parallel

mkdir forster
split -d -l 6 the_file_for_foster.fold forster/com
cd forster
ls *> fastaname1.txt
cat fastaname1.txt | parallel -j 5 --keep-order 'RNAforester -f {} -r > forester{}.S'
cat *.S > forester_S
rm *.S
cat fastaname1.txt | parallel -j 5 --keep-order 'RNAforester -f {} -d > forester{}.d'
cat forester*.d > forester_d
rm *.d

cd ..
mkdir result
mv forster/forester_S result
mv forster/forester_d result
rm -rf forster
cd result

grep "global optimal score: " forester_S | awk '{print $4}' > similarity_score.txt
grep -A 1 "global optimal score" forester_S | sed -n '/global optimal score/{n;p}' > relative_score.txt
grep -A 1 "global optimal score" forester_d | sed -n '/global optimal score/{n;p}' >site_name.txt
sed -i 's/[[:space:]].*//' site_name.txt
sed -i 's/.$//' site_name.txt
grep "global optimal score: " forester_d | awk '{print $4}' > forester_distance.txt

mv ../RNAstructure_fastam/MFE_structure_m6A.txt MFE_structure_m6A
mv ../RNAstructure_fasta/MFE_structure.txt MFE_structure
rmdir ../RNAstructure_fastam
rmdir ../RNAstructure_fasta
awk '{print substr($0, 11)}' MFE_structure_m6A > MFE_structure_m6A_R.txt
awk '{print substr($0, 11)}' MFE_structure > MFE_structure_R.txt
cd ..


#######file rename
directory_path="CTm"

cd "$directory_path" || exit

for file in *; do
    if [[ -f $file ]]; then  
        new_name=$(head -n 1 "$file" | sed -n 's/.*m6A\(.*\)/CTmm6A\1/p')  
        if [ -n "$new_name" ]; then  
            mv "$file" "$new_name"  
        fi
    fi
done
cd ..

directory_path="CTa"

cd "$directory_path" || exit

for file in *; do
    if [[ -f $file ]]; then  
        new_name=$(head -n 1 "$file" | sed -n 's/.*m6A\(.*\)/CTam6A\1/p') 
        if [ -n "$new_name" ]; then 
            mv "$file" "$new_name" 
        fi
    fi
done



