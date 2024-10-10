# -*- encoding: utf-8 -*-
'''
@File    :   compute_Fox60_bins2.py
@Time    :   2024/9/21 08:51:14
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

def dna_bins(coorddata):
    
    binwidth = 10
    dnal=100
    dnar=200
    boxl = 0
    boxr = 300
    coorddata = np.asfarray(coorddata)
    coorddata[:, 0] = coorddata[:, 0]+(coorddata[:, 0]//30)*(-30)
    coorddata[:, 1] = coorddata[:, 1]+(coorddata[:, 1]//30)*(-30)
    coorddata[:, 2] = coorddata[:, 2]+(coorddata[:, 2]//300)*(-300)
    coorddata = np.asfarray(coorddata).T

    # print(coorddata.shape)
    coorddata_z = coorddata[2]
    sections = np.arange(boxl, (boxr+0.1), binwidth)
    # sections[0] -= 5
    # sections[-1] += 5
    binsdividedata = pd.cut(coorddata_z, sections)
    binscounts = pd.value_counts(
        binsdividedata, sort=False).to_frame().reset_index()
    binscounts.columns = ['intervals', 'counts']
    # print(binscounts['intervals'])
    dnabins = np.array(binscounts['counts'])
    # print(dnabins)

    return dnabins


if __name__ == "__main__":
    avgn1 = 60*200
    fnames = locals()
    totalframe = 10001 # last 1us
    timeinterval = 1 # 0.1ns
    # max_filtration_value=4.0
    results =[]
    for index, inpfil in enumerate(sys.argv[1:]):
        inpdata1process = []
        # inpdata2process = []
        totalrownum1 = sum_count(inpfil)
        with open(inpfil, 'r') as f_handle:
            print(int(totalrownum1/avgn1))
            for nsampl in np.arange(int(totalrownum1/avgn1)):
                # for nsampl in np.arange(3446):
                onedata = read_fileblock(f_handle, avgn1)
                # countsampl += 1
                onedata = np.array(onedata)[:,1:]
                # onedata = np.asfarray(onedata)
                inpdata1process.append(dna_bins(onedata))
        inpdata1process = np.asfarray(inpdata1process)
        results.append(np.average(inpdata1process, axis=0))
        results.append(np.std(inpdata1process, ddof=1, axis=0))
    results = np.array(results)
    with open("/public/home/Skleac_wj_14/PDB/FOXA1_HP1a/analysis/0.82ePDle70/data_Fox60_bins2.txt", 'a') as savef:
        savef.write("# Fox60 in 30 bins \n")
        for line in results:
            savef.write(" ".join(str('{:.2f}'.format(value))
                        for value in line) + "\n")
    print("done\n")
    
