#!/bin/bash

RUNROOT=`pwd`
#Ingro="$RUNROOT/hp1a_idr.29567.pdb.sb/hp1a_idr.29567.pdb.gro"
#Ingro="$RUNROOT/foxa1_idr48aa.29566.pdb.sb/foxa1_idr48aa.29566.pdb.gro"
#Ingro="$RUNROOT/foxa1_idr.29565.pdb.sb/foxa1_idr.29565.pdb.gro"
Ingro="$RUNROOT/Fox60.29593.pdb.gro"
Ougro="Fox60"   #"foxa1" "foxa1_48aa"    #"hp1a"
len=60  #124  48
Mtype="FOX"  #"HPA"
Num=200 #100  200

#cat $Ingro | awk -v var=$len 'NR>=3 && NR<=(var+2)' > gro_tmp1
#cat gro_tmp1 | awk -v var=$Mtype '{printf "    1%5s%5s%5d %7.3f %7.3f %7.3f\n", var,$2,$4,$5,$6,$7}' > gro_tmp2
#head -2 $Ingro >  gro_tmp
#cat gro_tmp gro_tmp2 > ${Ougro}_seed.gro
#tail -1 $Ingro >> ${Ougro}_seed.gro
#
#
#rm gro_tmp1 gro_tmp2 gro_tmp

#gmx insert-molecules -ci ${Ougro}_seed.gro -nmol $Num -box 80 80 80 -o ${Ougro}_${Num}_box80.gro
#gmx editconf -f ${Ougro}_${Num}_box80.gro -o ${Ougro}_${Num}_box80_center.gro -center 0 0 0
gmx insert-molecules -ci ${Ougro}_seed.gro -nmol $Num -box 30 30 300 -o ${Ougro}_${Num}_box1.gro

#rm ${Ougro}_${Num}_box80.gro
#gmx editconf -f pld_CA_seed.gro -o pld_CA_1_box40_center.gro -c -d 20 -bt cubic

