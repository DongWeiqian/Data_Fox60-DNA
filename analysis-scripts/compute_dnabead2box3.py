# -*- encoding: utf-8 -*-
'''
@File    :   compute_dnabead2box3.py
@Time    :   2024/9/21 08:42:16
@Author  :   Dong Weiqian
@Version :   1.0
@Contact :   dwqian@ciac.ac.com
'''
### now start 


# import gudhi
import matplotlib
# Solve ssh GUI don't display the figure by python
matplotlib.use('tkagg') # Must be before importing matplotlib.pyplot or pylab!
import numpy as np
import matplotlib.pylab as plt
import sys
from typing import List
import copy
import pandas as pd
import math
from scipy.spatial import distance

def wc_count(file_name):
    import subprocess
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])

def sum_count(file_name):
    return sum(1 for _ in open(file_name))

def read_fileblock(filehandle, nline):
    read_datablock = []
    i = 0
    while filehandle:
        line = filehandle.readline()
        if not line:
            break
        line = ''.join(line)
        line = line.strip("\n")
        line = line.split()
        # print(line[0])
        if not line[0].startswith("#"):
            read_datablock.append(line)
            i += 1
            # if len(line) <= 3:
            #     print("error", read_datablock[0], line)
            #     break
            if i == nline:
                #print(filehandle.tell())
                break 
    return read_datablock


def dna_box(coorddata, highcenter):
    # Unit lenth=nm
    highcenterx = highcenter[0]
    highcentery = highcenter[1]
    highcenterz = highcenter[2]
    totalbeads = 980
    pchainlen = 60
    coorddata = np.asfarray(coorddata)
    # coorddata[:, 0] = coorddata[:, 0]+(coorddata[:, 0]//30)*(-30)
    # coorddata[:, 1] = coorddata[:, 1]+(coorddata[:, 1]//30)*(-30)
    # coorddata[:, 2] = coorddata[:, 2]+(coorddata[:, 2]//300)*(-300)
    d_highx = coorddata[:, 0]-highcenterx
    coorddata[:, 0] = d_highx+np.array(d_highx/15).astype(int)*(-30)
    d_highy = coorddata[:, 1]-highcentery
    coorddata[:, 1] = d_highy+np.array(d_highy/15).astype(int)*(-30)
    d_highz = coorddata[:, 2]-highcenterz
    coorddata[:, 2] = d_highz +np.array(d_highz/150).astype(int)*(-300)  # mention
    # if coorddata.shape[0] == 2 :
    #     print(coorddata[:,0])
    width = 5
    coorddata1_df = pd.DataFrame(coorddata)
    # print(comcoorddata_df.shape)
    # highcenterindex = comcoorddata1_df[(comcoorddata1_df[0] >= high_xf) & (comcoorddata1_df[0] <= high_xb) & (comcoorddata1_df[1] >= high_yd) & (comcoorddata1_df[1] <= high_yu) & (comcoorddata1_df[2] >= high_zleft) & (comcoorddata1_df[2] <= high_zright)].index.tolist()
    highcenterdata = coorddata1_df[(coorddata1_df[0] >= -width) & (coorddata1_df[0] <= width) & (coorddata1_df[1] >= -width) & (
        coorddata1_df[1] <= width) & (coorddata1_df[2] >= -width) & (coorddata1_df[2] <= width)]
    # print(np.array(highcenterindex, dtype=int).shape)
    num_dna_cbox = highcenterdata.shape[0]

    return num_dna_cbox


def compute_dnacenter(dnadata):
    coorddata = np.array(dnadata)
    binwidth = 10
    dnal = 100
    dnar = 200
    if len(coorddata)>0: 
        coorddata = np.asfarray(coorddata)
        # print(coorddata[0])
        # coorddata[:, 0] = coorddata[:, 0]+(coorddata[:, 0]//30)*(-30)
        # coorddata[:, 1] = coorddata[:, 1]+(coorddata[:, 1]//30)*(-30)
        # coorddata[:, 2] = coorddata[:, 2]+(coorddata[:, 2]//300)*(-300)
        bins = np.arange((dnal-10), (dnar+10.1), binwidth)
        df = pd.DataFrame(data=coorddata, columns=['x','y','z'])
        df['bin'] = pd.cut(df['z'], bins=bins)

        # Calculate bin counts
        bin_counts = df['bin'].value_counts()

        # Find the bin with the maximum count
        max_bin = bin_counts.idxmax()
        # find_bin = pd.cut([highcenterz], bins=bins)[0]

        # Output all data in the bin with the maximum count
        # find_bin_data = df[df['bin'] == find_bin]
        # if find_bin_data.shape[0] >= 10:
        #     max_bin_data = find_bin_data
        # else:
        #     max_bin_data = df[df['bin'] == max_bin]
        max_bin_data = df[df['bin'] == max_bin]
        coorddata_center = np.array([max_bin_data['x'].mean(),max_bin_data['y'].mean(),max_bin_data['z'].mean()])
        # pcenter = np.array([highcenterx,highcentery,highcenterz])
        # coorddata_center = np.average(coorddata, axis=0)
        # dist = np.linalg.norm((coorddata_center-pcenter))
        # print(coorddata_center)
        beadnum_box = dna_box(dnadata, coorddata_center)
        # print(beadnum_box)
    else: 
        beadnum_box = 1000
    

    return beadnum_box


if __name__ == "__main__":
    # avgn1 = 12000 
    # fnames = locals()
    totalframe = 10001 # last 1us
    timeinterval = 1 # 0.1ns
    box_x = 30
    box_y = 30
    box_z = 300
    # max_filtration_value=4.0
    results1 = []
    results2 = []
    
    # fig1, axes1 = plt.subplots(int(len(sys.argv)/2), 2,
    #                            figsize=(6, int(len(sys.argv)/2)*3))
    for index, inpfil in enumerate(sys.argv[1::2]):
        inpdata1process = []
        inpdata2process = []
        totalrownum = sum_count(inpfil)
        dnalen = 490
        interval = 6
        highcenterxyz = np.array(np.loadtxt(sys.argv[(index*2+2)], comments="#"))[:,:]
        # print(sys.argv[(index*2+2)])
        with open(inpfil, 'r') as f_handle:
            # print(np.arange(int(totalrownum/avgn1)))
            print(int(totalrownum/dnalen))
            for nsampl in np.arange(int(totalrownum/dnalen)):
            # for nsampl in np.arange(3446):
                onedata = read_fileblock(f_handle, dnalen)
                onedata = np.array(onedata)
                onedata = np.asfarray(onedata)
                onedata2 = [[i for i in j] for j in onedata]
                highcenterx = highcenterxyz[nsampl, 0]
                highcentery = highcenterxyz[nsampl, 1]
                highcenterz = highcenterxyz[nsampl, 2]
                if highcenterx < 0:
                    highcenterx = highcenterx + box_x
                elif highcenterx > box_x:
                    highcenterx = highcenterx - box_x
                if highcentery < 0:
                    highcentery = highcentery + box_y
                elif highcentery > box_y:
                    highcentery = highcentery - box_y
                if highcenterz < 0:
                    highcenterz = highcenterz + box_z
                elif highcenterz > box_z:
                    highcenterz = highcenterz - box_z
                
                pdnum = dna_box(onedata, [highcenterx, highcentery, highcenterz])
                
                inpdata1process.append(pdnum)
                # print(onedata2[0])
                inpdata2process.append(compute_dnacenter(onedata2))
                # print(inpdata2process)
                
        inpdata1process = np.array(inpdata1process)
        inpdata2process = np.array(inpdata2process)
        # print(inpdata1process.shape)
        # results.append([np.average(inpdata1process[:, 0], axis=0), np.std(inpdata1process[:, 0], ddof=1, axis=0), np.average(
        #     inpdata1process[:, 1], axis=0), np.std(inpdata1process[:, 1], ddof=1, axis=0)])
        results1.append([np.average(inpdata1process, axis=0), np.std(inpdata1process, ddof=1, axis=0)])
        results2.append([np.average(inpdata2process, axis=0), np.std(inpdata2process, ddof=1, axis=0)])  
    results1 = np.asfarray(results1).T
    results2 = np.asfarray(results2).T
    # print(results)
    # results = [0,1,2,3]
    with open("/public/home/Skleac_wj_14/PDB/FOXA1_HP1a/analysis/0.82ePDle70/dnabead2box.txt", 'a') as savef:
        savef.write("# dnabead 490beads number in protein center box in dna max bead bin center box \n")
        savef.write(" ".join(str('{:.2f}'.format(value)) for value in results1[0]) + "\n")
        savef.write(" ".join(str('{:.2f}'.format(value)) for value in results1[1]) + "\n")
        savef.write(" ".join(str('{:.2f}'.format(value)) for value in results2[0]) + "\n")
        savef.write(" ".join(str('{:.2f}'.format(value)) for value in results2[1]) + "\n")
    print("done\n")
    
