# -*- encoding: utf-8 -*-
'''
@File    :   compute_proteinsrho.py
@Time    :   2024/02/24 18:11:07
@Author  :   Dong Weiqian
@Version :   1.6
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
# Compute the highest density, lowest density and density difference of proteins along z axis
def computerho_highlowdiffer(coorddata, highcenterx, highcentery, highcenterz, box_z, nsample):
    # Unit lenth=nm
    totalbeads = 12000
    binwidth = 10
    coorddata = np.asfarray(coorddata).T
    # if coorddata.shape[0] == 2 :
    #     print(coorddata[:,0])
    coorddata[0] = coorddata[0]+(coorddata[0]//30)*(-30)
    coorddata[1] = coorddata[1]+(coorddata[1]//30)*(-30)
    coorddata[2] = coorddata[2]+(coorddata[2]//300)*(-300)
    d_highx = coorddata[0]-highcenterx
    coorddata[0] = d_highx+np.array(d_highx/15).astype(int)*(-30)
    d_highy = coorddata[1]-highcentery
    coorddata[1] = d_highy+np.array(d_highy/15).astype(int)*(-30)
    coorddata_z = coorddata[2]
    highwidth = 10
    high_zleft = highcenterz - highwidth/2
    high_zright = highcenterz + highwidth/2
    transitionwidth = 50
    dnalend = 100
    dnarend = 170
    low_zll = 10
    low_zlr = dnalend - transitionwidth
    low_zrl = dnarend + transitionwidth
    low_zrr = box_z - 10

    sections = [10, 50, 220, 290]
    group_names = ["L1", "H", "L2"]
    binsdividedata = pd.cut(coorddata_z, sections, labels=group_names)
    binscounts = pd.value_counts(binsdividedata).sort_index()
    value_counts = binscounts.values.tolist()
    rho_low = (value_counts[0]+value_counts[2]) / (totalbeads*(50-10+290-220)/binwidth)
    # if (high_zleft > low_zlr) and (high_zright < low_zrl):
    #     sections = [low_zll, low_zlr, high_zleft, high_zright, low_zrl, low_zrr]
    #     # if (high_zleft < low_zlr) or (high_zright > low_zrl):
    #     #     print("Fatal Error: Group arange is", sections)
    #     group_names = ["L1", "M1", "H", "M2", "L2"]
    #     binsdividedata = pd.cut(coorddata_z, sections, labels=group_names)
    #     binscounts = pd.value_counts(binsdividedata).sort_index()
    #     # print(binscounts)
    #     value_counts = binscounts.values.tolist()
    #     # rho_high = value_counts[2]/(totalbeads*highwidth/binwidth)
    #     rho_low = (value_counts[0]+value_counts[4])/(totalbeads*(low_zrr-low_zrl+low_zlr-low_zll)/binwidth)
    #     # rho_differ = rho_high - rho_low
    # elif (high_zleft < low_zlr) and (high_zright < low_zrl):
    #     sections = [high_zleft, high_zright, low_zrl, low_zrr]
    #     group_names = ["H", "M1", "L1"]
    #     binsdividedata = pd.cut(coorddata_z, sections, labels=group_names)
    #     binscounts = pd.value_counts(binsdividedata).sort_index()
    #     # print(binscounts)
    #     value_counts = binscounts.values.tolist()
    #     # rho_high = value_counts[0]/(totalbeads*highwidth/binwidth)
    #     rho_low = value_counts[2] / (totalbeads*(low_zrr-low_zrl)/binwidth)
    #     # rho_differ = rho_high - rho_low
    # elif (high_zleft > low_zlr) and (high_zright > low_zrl):
    #     sections = [low_zll, low_zlr, high_zleft, high_zright]
    #     group_names = ["L1", "M1", "H"]
    #     binsdividedata = pd.cut(coorddata_z, sections, labels=group_names)
    #     binscounts = pd.value_counts(binsdividedata).sort_index()
    #     # print(binscounts)
    #     value_counts = binscounts.values.tolist()
    #     # rho_high = value_counts[2]/(totalbeads*highwidth/binwidth)
    #     rho_low = value_counts[0] / (totalbeads*(low_zlr-low_zll)/binwidth)
    #     # rho_differ = rho_high - rho_low
    # else:
    #     print(highcenterz,high_zleft,high_zright,nsample)
    # print(rho_high, rho_low, rho_differ, highcenter)
    d_highz = coorddata[2]-highcenterz
    coorddata[2] = d_highz+np.array(d_highz/150).astype(int)*(-300) #mention
    width =5
    high_xf = highcenterx-5
    high_xb = highcenterx+5
    high_yd = highcentery-5
    high_yu = highcentery+5
    coorddata_df = pd.DataFrame(coorddata.T)
    # print(coorddata_df.shape)
    # highcenterdata = coorddata_df[(coorddata_df[0] >= high_xf) & (coorddata_df[0] <= high_xb) & (coorddata_df[1] >= high_yd) & (coorddata_df[1] <= high_yu) & (coorddata_df[2] >= high_zleft) & (coorddata_df[2] <= high_zright)]
    # highcenterdata = coorddata_df[(coorddata_df[0] >= -width) & (coorddata_df[0] <= width) & (coorddata_df[1] >= -width) & (
    #     coorddata_df[1] <= width) & (coorddata_df[2] >= -(width+5)) & (coorddata_df[2] <= (width+5))]
    highcenterdata = coorddata_df[(coorddata_df[0] >= -width) & (coorddata_df[0] <= width) & (coorddata_df[1] >= -width) & (
        coorddata_df[1] <= width) & (coorddata_df[2] >= -(width)) & (coorddata_df[2] <= (width))]
    # print(highcenterx,highcentery,highcenterdata.shape)
    rho_high = highcenterdata.shape[0]/(totalbeads*highwidth/binwidth) * 1947.987 # mg/ml
    rho_low = rho_low/9 * 1947.987
    rho_differ = rho_high - rho_low
    return rho_high, rho_low, rho_differ

def excutemultifile(index, inpfil, parameter1data, timeinterval):
    inpdata1process = []
    # inpdata2process = []
    totalrownum = sum_count(inpfil)
    avgn1 = 60*200
    # avgn2 = 1000 # 100ns
    avgn2 = timeinterval
    box_x = 30
    box_y = 30
    box_z = 300
    tmp1data = []
    # tmp2data = []
    parameter1avg = []
    with open(inpfil, 'r') as f_handle:
        # print(np.arange(int(totalrownum/avgn1)))
        countsampl = 0
        print(int(totalrownum/avgn1))
        for nsampl in np.arange(int(totalrownum/avgn1)):
        # for nsampl in np.arange(3446):
            onedata = read_fileblock(f_handle, avgn1)
            countsampl += 1
            onedata = np.array(onedata)[:, 1:]
            onedata = np.asfarray(onedata)
            highcenterx = parameter1data[nsampl,0]
            highcentery = parameter1data[nsampl,1]
            highcenterz = parameter1data[nsampl,2]
            
            if np.isnan(highcenterx):
                continue
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
                # 
            # print(highcenterx,highcentery,highcenterz, nsampl)
            # pre_centerxyz = parameter1data[nsampl,0:]
            rho_highlowdiffer = computerho_highlowdiffer(onedata, highcenterx, highcentery, highcenterz, box_z, nsampl)
            inpdata1process.append(rho_highlowdiffer)
            
    with open("rho_tgridhighlowdiffer3", 'w') as savef:
        savef.write("# rho_high_gridxyz  rho_low  rho_differ *1947.987mg/ml\n")
        if inpdata1process: 
            for line in inpdata1process:
                savef.write(" ".join(str('{:.5f}'.format(value)) for value in line) + "\n")
    return index, inpfil, inpdata1process, parameter1avg

if __name__ == "__main__":
    avgn1 = 12000 
    fnames = locals()
    totalframe = 10001 # last 1us
    timeinterval = 1 # 0.1ns
    # max_filtration_value=4.0
    results =[]
    parameter1data = np.array(np.loadtxt(sys.argv[-1], comments="#"))[:,:]
    
    for index, inpfil in enumerate([sys.argv[1]]):
       results.append(excutemultifile(index, inpfil, parameter1data, timeinterval)) 
    print("done\n")
    