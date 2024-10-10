# -*- encoding: utf-8 -*-
'''
@File    :   compute_dnafc2.py
@Time    :   2024/9/21 08:45:29
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
    dnafoldcontanctnum = 0
    rd = 1.5
    for ni in np.arange((np.array(dnadata).shape[0]-interval-1)):
        nj = int(ni+interval+1)
        nidata = alldist[ni, nj:]
        nidata[nidata <= rd ] = 1
        nidata[nidata > rd] =0
        dnafoldcontanctnum = dnafoldcontanctnum + np.sum(nidata)

    return dnafoldcontanctnum


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
                dnafoldcontact = compute_dnafoldcontact(onedata,interval)
                inpdata1process.append(dnafoldcontact)
        inpdata1process = np.array(inpdata1process)
        results.append([np.average(inpdata1process, axis=0),np.std(inpdata1process, ddof=1, axis=0)])
    results = np.asfarray(results).T
    
    with open("/public/home/Skleac_wj_14/PDB/FOXA1_HP1a/analysis/0.82ePDle70/dnafc2all2.txt", 'a') as savef:
        savef.write("# dnafoldcontact rd=1.5 490beads \n")
        # savef.write(str('{:.3f}'.format(np.average(inpdata1process, axis=0))) + "  ")
        # savef.write(str('{:.3f}'.format(np.std(inpdata1process, ddof=1))) + "\n")
        savef.write(" ".join(str('{:.2f}'.format(value)) for value in results[0]) + "\n")
        savef.write(" ".join(str('{:.2f}'.format(value)) for value in results[1]) + "\n")
    
    print("done\n")
    
