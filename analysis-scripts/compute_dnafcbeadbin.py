# -*- encoding: utf-8 -*-
'''
@File    :   compute_dnafcbeadbin.py
@Time    :   2024/9/21 08:48:25
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
            if i == nline:
                #print(filehandle.tell())
                break 
    return read_datablock

def compute_dnafoldcontact(dnadata, interval):
    # pos_ele = [1,1,1,1,1,1,1,1,1,0.25,1,0.25,1]
    alldist = distance.cdist(dnadata, dnadata)
    alldist = np.array(alldist)
    beadnum = 0
    contactnum = []
    rd = 1.5
    internum = 7
    rchange = np.array(dnadata).shape[0]-interval-1
    for ni in np.arange(490):
        nj = int(ni+internum)
        if (ni >= internum) and (ni <= (489-internum)):
            nidata = np.append(alldist[ni, :(ni-internum)], alldist[ni, nj:])
            nidata = np.array(nidata)
        elif ni < internum:
            nidata = alldist[ni, nj:]
        else:
            nidata = alldist[ni, :(ni-internum)]
        nidata[nidata <= rd ] = 1
        nidata[nidata > rd] =0
        # print(nidata)
        if np.sum(nidata) >=1 :
            contactnum.append(ni)
    coorddata = np.array(dnadata)[contactnum]
    binwidth = 10
    dnal = 100
    dnar = 200
    coorddata = np.asfarray(coorddata)
    coorddata[:, 0] = coorddata[:, 0]+(coorddata[:, 0]//30)*(-30)
    coorddata[:, 1] = coorddata[:, 1]+(coorddata[:, 1]//30)*(-30)
    coorddata[:, 2] = coorddata[:, 2]+(coorddata[:, 2]//300)*(-300)
    coorddata = np.asfarray(coorddata).T

    # print(coorddata.shape)
    coorddata_z = coorddata[2]
    sections = np.arange((dnal-10), (dnar+10.1), binwidth)
    # sections[0] -= 5
    # sections[-1] += 5
    binsdividedata = pd.cut(coorddata_z, sections)
    binscounts = pd.value_counts(
        binsdividedata, sort=False).to_frame().reset_index()
    binscounts.columns = ['intervals', 'counts']
    # print(binscounts['intervals'])
    dnabins = np.array(binscounts['counts'])
    return dnabins


if __name__ == "__main__":
    # avgn1 = 12000 
    # fnames = locals()
    totalframe = 10001 # last 1us
    timeinterval = 1 # 0.1ns
    # max_filtration_value=4.0
    results =[]
    for index, inpfil in enumerate(sys.argv[1:]):
        inpdata1process = []
        # inpdata2process = []
        totalrownum = sum_count(inpfil)
        dnalen = 490
        interval = 6
        with open(inpfil, 'r') as f_handle:
            # print(np.arange(int(totalrownum/avgn1)))
            print(int(totalrownum/dnalen))
            for nsampl in np.arange(int(totalrownum/dnalen)):
            # for nsampl in np.arange(3446):
                onedata = read_fileblock(f_handle, dnalen)
                onedata = np.array(onedata)
                onedata = np.asfarray(onedata)
                beadnum = compute_dnafoldcontact(onedata,interval)
                inpdata1process.append(beadnum)
        inpdata1process = np.array(inpdata1process)
        results.append(np.average(inpdata1process, axis=0))
        results.append(np.std(inpdata1process, ddof=1, axis=0))
    results = np.asfarray(results)
    with open("/public/home/Skleac_wj_14/PDB/FOXA1_HP1a/analysis/0.82ePDle70/dnafcbeadbin.txt", 'a') as savef:
        savef.write("# dnafoldcontact bead number in 12 bins rd=1.5 490beads \n")
        for line in results:
            savef.write(" ".join(str('{:.2f}'.format(value))
                        for value in line) + "\n")
        # savef.write(" ".join(str('{:.2f}'.format(value)) for value in results[0]) + "\n")
        # savef.write(" ".join(str('{:.2f}'.format(value)) for value in results[1]) + "\n")
    print("done\n")
    
