#!/bin/bash
# dwq 2024/9/3
#Nam="Fox60le70"
Job="box1le70Fox60"
Xbox=30

Inanaly="/public/home/Skleac_wj_14/PDB/FOXA1_HP1a/analysis/0.82ePDle70"

for i in 250 270 275 280 285 290 300 310
do
cd temp${i}
awk 'NF==4 && $1=="X" {print $2, $3, $4}' adj7interv6pc1_center.xyz > adj7pc1_center.xyz
awk 'NF==4 && $1=="X" {print $2, $3, $4}' proteins_com.xyz > p200_com.xyz

bsub -q serial -n 1 -J ${i}rho python $Inanaly/compute_proteinsrho_xyzgrid2.py last1us200proteins.xyz adj7pc1_center.xyz

cd ..
printf "${i} \n"
done

fname="arrsdna.xyz"
pyexe="compute_dnafcbeadbin.py"
python $Inanaly/$pyexe temp250/$fname temp270/$fname temp275/$fname 

pyexe="compute_dnafc2.py"
python $Inanaly/$pyexe temp250/$fname temp270/$fname temp275/$fname 

fname="last1us200proteins.xyz"
pyexe="compute_Fox60_bins2.py"
python $Inanaly/$pyexe temp250/$fname temp270/$fname temp275/$fname 

fname="protednadist1.data"
pyexe="print_dnasegprotecont.py"
python $Inanaly/$pyexe temp250/$fname temp270/$fname temp275/$fname 

pyexe="compute_dnabead2box3.py"
fname="arrsdna.xyz"
fname2="adj7pc1_center.xyz"
python $Inanaly/$pyexe temp250/$fname temp250/$fname2 temp270/$fname temp270/$fname2 temp275/$fname temp275/$fname2


