#!/usr/bin/env python3

from hcalConfig import *
import math,sys

pmts=readPMTList(read_summary=True)
norms=readNorms()
max_npe=7805.99


def printPMTVal(pmt,sv,run):
    if not sv:
        return
    norm=norms[sv.row][sv.col]
    norm_npe=sv.npe/norm
    qe=norm_npe/max_npe
    print('%7s %4d %2d %2d %7.1f %7.3f %7.1f %7.3f'%(pmt,run,sv.row,sv.col,sv.npe,norm,norm_npe,qe))

if len(sys.argv) == 2:
    #run = next((r for r in runs if sv.row == row and sv.col == col and sv.led==32),None)
    for run,l in runs.items():
        pmt=next((ll for ll in l if ll[2] == sys.argv[1]),None)
        if pmt != None:
            printPMTVal(sys.argv[1],pmt[3],run)
else:
    for pmt,l in pmts.items():
        printPMTVal(pmt,l[3],l[2])

