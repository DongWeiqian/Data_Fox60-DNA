# -*- encoding: utf-8 -*-
'''
@File    :   print_dnaprotecont.py
@Time    :   2024/08/07 15:17:43
@Author  :   Dong Weiqian
@Version :   1.0
@Contact :   dwqian@ciac.ac.com
'''
### now start 

import sys
import matplotlib.pylab as plt
import numpy as np
from adjustText import adjust_text
import matplotlib
# Solve ssh GUI don't display the figure by python
matplotlib.use('tkagg')  # Must be before importing matplotlib.pyplot or pylab!
# Count big file total lines
from scipy.optimize import curve_fit

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
                # print(filehandle.tell())
                break
    return read_datablock

if __name__ == "__main__":
    avgn1 = 10000
    data_mean1 = []
    data_mean2 = []
    #print(-17+int(-17/15)*(-30),17+int(17/15)*(-30),-12+int(-12/15)*(-30),12+int(12/15)*(-30))
    # fnames = locals()
    for index, inpfil in enumerate(sys.argv[1:]):
        # fnames["inpdata"+str(index+1)] = []
        with open(inpfil, 'r') as f_handle:
            # for nsampl in np.arange(int(wc_count(inpfil)/avgn1)):
            for nsampl in np.arange(1):
                # print(wc_count(inpfil))
                onedata = read_fileblock(f_handle, int(wc_count(inpfil)))
                onedata1 = np.asfarray(onedata)[:, 1:]  # DNA 10-bead seg 490beads
                data_mean1.append(np.average(onedata1, axis=0))
                data_mean1.append(np.std(onedata1, ddof=1, axis=0))
    data_mean1 = np.asfarray(data_mean1)
    # data_mean2 = np.asfarray(data_mean2).T
    print(np.sum(data_mean1[0]), np.sum(data_mean1[2]), np.sum(data_mean1[4]))
    with open("/public/home/Skleac_wj_14/PDB/FOXA1_HP1a/analysis/0.82ePDle70/dnaseg_Fox60_cont1.txt", 'a') as savef:
        savef.write("# number of protein beads contacting with DNA 10-beads seg 490 beads rd=1.5 nm \n")
        
        for line in data_mean1:
            savef.write(" ".join(str('{:.2f}'.format(value))
                        for value in line) + "\n")
        # savef.write(" ".join(str('{:.2f}'.format(value)) for value in data_mean1[0]) + "\n")
        # savef.write(" ".join(str('{:.2f}'.format(value)) for value in data_mean1[1]) + "\n")
