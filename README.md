# Fox60-DNA LLPS MD Guide
By Weiqian Dong
----------------------

1. We have shown the data of figures used in our research article (Regulatory Mechanisms of Nanoscale FOXA1 and DNA Condensates Under DNA Tension), like Fig1.Disorder_score_IUPred3_1.txt.
2. We have prepared the intial coarse-grained Fox60-DAN model for LLPS MD in MOFF-DNA_CGmodel folder. If you want to generate your coarse-grained model or know the detail information of generating MOFF-DNA model, you can refer Latham, A.P. and B. Zhang, On the stability and layered organization of protein-DNA condensates. Biophysical Journal, 2022. 121(9): p. 1727-1737 (https://doi.org/10.1016/j.bpj.2022.03.028).
3. In analysis-scripts folder, we shown the codes of analyzing our MD trajectory data.
   -- File interval6proteincluster-adj7.dat is used for calculating N_clus and finding the center of protein condensates by PLUMED 2.8.
   -- File compute_Fox60_bins2.py and compute_dnafcbeadbin.py are used for calculating $B_{pro\{_}z}$ and $B_{contdna\\_z}$, respectively.
   -- File compute_proteinsrho_xyzgrid2.py is used for calculating $\Delta\rho$.
   -- File compute_dnafc2.py is used for calculating $N_{dna\_intra}$.
   -- File compute_dnabead2box3.py is used for calculating $\rho_{dp}$ and $\rho_{dd}$.
   -- File print_dnasegprotecont.py is used for printing out $N_{PD\_inter}$.
   -- File protednadist1.dat is used for generating the inter-chain contact number between DNA and protein beads by PLUMED 2.8.
   -- File sdna_xyz.dat is used for generating the coordinate of a single $DNA_{490}$ by PLUMED 2.8.
   -- File runanaly1.sh is used for coverting the MD trajectory files (.xtc files) to readable matrix files and preparing the input files for computing data analysis.
   -- File compute1.sh is sh command used for computing data analysis.

------------------------
Email dwqian@ciac.ac.com with any questions
