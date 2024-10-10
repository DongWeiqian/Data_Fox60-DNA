#!/bin/bash
# dwq 2024/4/11
#Nam="Fox60le70"
Job="box1le70Fox60"
Xbox=30

Infolder="/public/home/Skleac_wj_14/PDB/FOXA1_HP1a/md/Fox60_DNA/box1le70"
Inanaly="/public/home/Skleac_wj_14/PDB/FOXA1_HP1a/analysis/0.82ePDle70"


for i in 250 270 275 280 285 290 300 310
do
rm -fr temp${i}
mkdir temp${i}
cd temp${i}
#rm -fr ${i}mM
#mkdir ${i}mM
#cd ${i}mM
#rm -fr 0e${i}
#mkdir 0e${i}
#cd 0e${i}
echo 2 | gmx trjconv -f $Infolder/temp${i}/step3_50mM_4fs_1.xtc -s $Infolder/temp${i}/step3_50mM_4fs_1.tpr -n $Inanaly/dna.ndx -b 4000000 -o tmpprote.gro -pbc whole > /dev/null 2>&1
awk '$1 ~ /FOX/ {print $1, $(NF-2), $(NF-1), $NF}' tmpprote.gro > last1us200proteins.xyz
rm tmpprote.gro
gmx trjconv -f $Infolder/temp${i}/step3_50mM_4fs_1.xtc -b 4000000 -o last1us.xtc

bsub -q serial -n 4 -J ${Job}${i}prote plumed driver --plumed $Inanaly/interval6proteincluster-adj7.dat  --ixtc last1us.xtc
bsub -q serial -n 1 -J ${i} plumed driver --plumed $Inanaly/protednadist1.dat --ixtc last1us.xtc
bsub -q serial -n 1 -J prote_com${i} plumed driver --plumed $Inanaly/prote_com.dat --ixtc last1us.xtc --mc $Inanaly/mass_charge_file

plumed driver --plumed $Inanaly/sdna_xyz.dat --ixtc last1us.xtc > /dev/null 2>&1
awk 'NF==4 && $1=="X" {print $2, $3, $4}' tmpsdna.xyz > arrsdna.xyz
rm tmpsdna.xyz
cd ..
printf "${i} \n"
done



