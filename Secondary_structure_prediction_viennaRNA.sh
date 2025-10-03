#!/bin/bash
#!/bin/bash
#SBATCH --job-name=Fold_Species
#SBATCH --partition=cpu6348
#SBATCH -n 5
#SBATCH --qos cpu6348
#SBATCH --mail-type=ALL
#SBATCH --output=%j.out
#SBATCH --error=%j.err

ml load rnastructure
ml load viennarna
ml load parallel
ml load r

export DATAPATH=/gpfs/spack/opt/linux-rocky8-icelake/gcc-8.5.0/rnastructure-6.4-yyu25welaovp5esuz6jafgl5fsi3zvct/data_tables

#######################The first part is for unmodified sequence 
input_file="pre_Species_sequence.fa" ###This is an example input
output_dir="RNAstructure_fasta"

mkdir -p $output_dir

awk -v RS=">" -v ORS="" '{if (NF) {print ">"$0 > "'$output_dir'/sequence_"++count".fasta"}}' $input_file

echo "FASTA file split into individual sequence files in the '$output_dir' directory."

cd RNAstructure_fasta

ls *.fasta > fastaname1.txt
cat fastaname1.txt | awk -F"." '{print $1}' > fastaname.txt
cat fastaname.txt | parallel -j 5 RNAfold --noPS  {}.fasta >  combined_seq 
cat fastaname.txt | parallel -j 5 rm {}.fasta
awk 'NR%3 == 0 {print $2}' combined_seq > MFE.txt
awk 'NR%3 == 1 {print $1}' combined_seq > site_name_MFE.txt
paste MFE.txt site_name_MFE.txt > MFE1.txt
sed 's/[>()]//g' MFE1.txt > MFE_structure_R.txt
rm MFE1.txt
rm fastaname1.txt
rm fastaname.txt
rm MFE.txt
rm site_name_MFE.txt
sed '0~3 s/ .*$//' combined_seq > combined_seq1
awk 'NR % 3 == 1 {$0 = $0"S"} {print}' combined_seq1 > RNAstructureS
rm combined_seq
awk '/^>/{if(x){close(x)} x=$0; sub("^>", "", x);}{print > (x".fold");}' combined_seq1

mkdir ../folda
ls *.fold > fastaname1.txt
cat fastaname1.txt| parallel -j 5 mv {} ../folda

mkdir ../combine
mv RNAstructureS ../combine

cd ../folda
for file in *; do mv "$file" "CTa$file"; done

cd ..

#########################################################
#######################The second part is for modified sequence 
input_file="pre_Species_sequence_modification.fa" ###This is an example input
output_dir="RNAstructure_fastam"

mkdir -p $output_dir

awk -v RS=">" -v ORS="" '{if (NF) {print ">"$0 > "'$output_dir'/sequence_"++count".fasta"}}' $input_file

echo "FASTA file split into individual sequence files in the '$output_dir' directory."

cd RNAstructure_fastam

ls *.fasta > fastaname1.txt
cat fastaname1.txt | awk -F"." '{print $1}' > fastaname.txt
cat fastaname.txt | parallel -j 5 RNAfold --noPS --modifications="I" {}.fasta >  combined_seq 
cat fastaname.txt | parallel -j 5 rm {}.fasta
awk 'NR%3 == 0 {print $2}' combined_seq > MFE.txt
awk 'NR%3 == 1 {print $1}' combined_seq > site_name_MFE.txt
paste MFE.txt site_name_MFE.txt > MFE1.txt
sed 's/[>()]//g' MFE1.txt > MFE_structure_m6A_R.txt
rm MFE1.txt
rm fastaname1.txt
rm fastaname.txt
rm MFE.txt
rm site_name_MFE.txt
sed '0~3 s/ .*$//' combined_seq > combined_seq1
sed '2~3s/I/A/g' combined_seq1 > combined_seq2
awk 'NR % 3 == 1 {$0 = $0"M"} {print}' combined_seq2  > RNAstructure_m6A

rm combined_seq
awk '/^>/{if(x){close(x)} x=$0; sub("^>", "", x);}{print > (x".fold");}' combined_seq2


mkdir ../foldm
ls *.fold > fastaname1.txt
cat fastaname1.txt| parallel -j 5 mv {} ../foldm

mv RNAstructure_m6A ../combine

cd ../foldm
for file in *; do mv "$file" "CTm$file"; done

cd ../combine
cat RNAstructureS RNAstructure_m6A > combined.fold

rm RNAstructureS
rm RNAstructure_m6A

###########################################################################################
####This is for RNAforester scoring 
awk '/^>/{if(x){close(x)} x=$0; sub("^>", "", x);}{print > x;}' combined.fold
rm *.fold
cat * > the_file_for_foster.fold
mv the_file_for_foster.fold ..
rm *

cd ..
rmdir combine


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

mv ../RNAstructure_fastam/MFE_structure_m6A_R.txt MFE_structure_m6A_R.txt
mv ../RNAstructure_fasta/MFE_structure_R.txt MFE_structure_R.txt 
rm -rf ../RNAstructure_fastam
rm -rf  ../RNAstructure_fasta
cd ..

mkdir CTa
mkdir CTm

cd folda
ls *.fold > fastaname1.txt
cat fastaname1.txt | awk -F"." '{print $1}' > fastaname.txt
cat fastaname.txt | parallel -j 5 dot2ct {}.fold {}.ct
mv *.ct ../CTa
rm fastaname1.txt
rm fastaname.txt

cd ../foldm
ls *.fold > fastaname1.txt
cat fastaname1.txt | awk -F"." '{print $1}' > fastaname.txt
cat fastaname.txt | parallel -j 5 dot2ct {}.fold {}.ct
mv *.ct ../CTm
rm fastaname1.txt
rm fastaname.txt

