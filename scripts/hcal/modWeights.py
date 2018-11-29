#!/usr/bin/env python3
use_col=3
scales = [1., 0.697, 1.]
#scales = [1., 1.0, 1.]

norm_to_row=None
norm_to_row=6 ## For consistency, normalize everything to row 6

configfile='configs/sub1col%d.cfg'%(use_col)
scale = scales[use_col-1]
printHighlight=True
errorWeightedAverage=1

from hcalConfig import *
import math

class StatList(object):
    __slots__=('vals','avg','std','min','max','n','rvals','errs','err')
    def __init__(self,vals=[],errs=[]):
        self.vals = list(vals)
        self.errs = list(errs)
        self.rvals = []
    def append(self,val,err):
        self.vals.append(val)
        self.errs.append(err)
    def Finalize(self):
        self.min =  1.e6
        self.max = -1.e6
        self.avg = 0.
        self.std = 0.
        self.err = 0.
        self.n = len(self.vals)
        if self.n <= 0:
            return
        elif self.n == 1:
            self.avg = vals[0]
            self.min = vals[0]
            self.max = vals[0]
            self.rvals.append(1.0)
            return
        for i in range(len(self.vals)):
            v = self.vals[i]
            ve = self.errs[i]*self.errs[i]
            if errorWeightedAverage:
                self.avg += v/ve
                self.err += 1./ve
            else:
                self.avg += v
            if self.min > v:
                self.min = v
            if self.max < v:
                self.max = v
        if errorWeightedAverage:
            self.avg /= self.err
            self.err = 1./math.sqrt(self.err)
        else:
            self.avg /= float(self.n)
        for v in self.vals:
            self.std  += (self.avg-v)*(self.avg-v)
            self.rvals.append(v/self.max)
        self.std = math.sqrt(self.std/float(self.n-1.))

class PMTNPE(object):
    __slots__=('pmt','npe', 'err','ratio','aratio','cal_npe')
    def __init__(self,pmt,npe,err):
        self.pmt = pmt
        self.npe = npe
        self.err = err
        self.ratio = 0.0
        self.aratio = 0.0
        self.cal_npe = 0.0
    def __lt__(self,other):
        return self.pmt < other.pmt

class RowCol(object):
    __slots__=('mod','row','col','pmts','max')
    @staticmethod
    def modnum(r,c):
        return (r-1)*12+(c-1)
    def __init__(self,row,col,pmts=[]):
        self.mod = self.modnum(row,col)
        self.row = row
        self.col = col
        self.pmts = list(pmts)
        self.max = None
    def Finalize(self):
        for p in self.pmts:
            if self.max == None or self.max.npe < p.npe:
                self.max = p
        if self.max == None:
            return

        ## Now that we established the maximum, take the ratios of all of them
        for p in self.pmts:
            p.ratio = p.npe/self.max.npe

def fmtVal(val,n=7,d=1,hi=False):
    ret = '  %'+'%d'%n+'.%d'%d+'f'
    if hi and printHighlight:
        ret += '*'
    else:
        ret += ' '
    return ret%(val)

runlist=readRunList(configfile)

## Parse the runlist
aMaxPMT=None
aMaxMod=None
runs=[]
modules=[]
first_run=True
nmodules=0
hdrmods1 = '  PMT: '
hdrmods2 = '-------'
for run,vals in runlist.items():
    summary = list(filter(lambda x: x.led == 32,readSummaryFile(run)))
    ## Loop through values and fill module dictionary accordingly
    for v in vals:
        ## Find module that corresponds to this
        row,col,pmt = v[0],v[1],v[2]
        idx = RowCol.modnum(row,col)
        mod = next((mod for mod in modules if mod.mod == idx),None)
        if mod == None:
            mod = RowCol(row,col)
            modules.append(mod)
        ## Now find the value in the summary list of this module
        idx += 32000
        sv = next((sv for sv in summary if sv.idx == idx))
        if sv != None:
            mod.pmts.append(PMTNPE(pmt,sv.npe,sv.err))
            if aMaxPMT is None or sv.npe > aMaxPMT.npe:
                aMaxPMT = PMTNPE(pmt,sv.npe,sv.err)
                aMaxModule = mod
    runs.append(run)

s_npe_mod=[]
s_npe_pmt=[]
s_rat_mod=[]
s_rat_pmt=[]
s_row_rat=[]
s_max_rat=[]


for im,m in enumerate(modules):
    m.pmts = sorted(m.pmts)
    m.Finalize()
    s_npe_mod.append(StatList())
    s_rat_mod.append(StatList())
    if im == 0:
        for p in m.pmts:
            hdrmods1 += '  %7s '%p.pmt
            hdrmods2 += '  ------- '
            s_npe_pmt.append(StatList())
            s_rat_pmt.append(StatList())
            s_row_rat.append(StatList())
            s_max_rat.append(StatList())

nruns=len(runs)

print(hdrmods1 + '    AVG       STDV ')
print(hdrmods2 + '   ------    ------')
for im,m in enumerate(modules):
    ln1 = '[%02d-%02d]'%(m.row,m.col)
    ln2 = '/MAXROW'
    ln3 = '/MAXALL'
    for ipmt,pmt in enumerate(m.pmts):
        #ln1 += '  %7s '%(pmt.pmt)
        is_max=False
        ln1 += fmtVal(pmt.npe,hi=(pmt.pmt==m.max.pmt))
        ln2 += fmtVal(pmt.ratio,d=2)
        pmt.aratio=pmt.npe/aMaxPMT.npe
        #pmt.aratio=pmt.npe/maxmaxnpe
        ln3 += fmtVal(pmt.aratio,d=3)
        s_npe_mod[im].append(pmt.npe,pmt.err)
        s_rat_mod[im].append(pmt.ratio,1.0)
        s_npe_pmt[ipmt].append(pmt.npe,1.0)
        s_rat_pmt[ipmt].append(pmt.ratio,1.0)
    ## Since we are done with module, we can set stats for that now
    s_npe_mod[im].Finalize()
    s_rat_mod[im].Finalize()
    ln1 += fmtVal(s_npe_mod[im].avg)
    ln2 += fmtVal(s_rat_mod[im].avg,d=2)
    ln1 += fmtVal(s_npe_mod[im].std)
    ln2 += fmtVal(s_rat_mod[im].std,d=2)

    print(ln1)
    print(ln2)
    print(ln3)
    print()

def printPMTStats(msg,lst,n=7,d=1):
    ## Now, compute the statistics for each PMT over all modules
    ln1 = 'AVG:   '
    ln2 = 'STD:   '
    ln3 = 'MIN:   '
    ln4 = 'MAX:   '
    for i in range(len(s_npe_pmt)):
        lst[i].Finalize();
        ln1 += fmtVal(lst[i].avg,n,d)
        ln2 += fmtVal(lst[i].std,n,d)
        ln3 += fmtVal(lst[i].min,n,d)
        ln4 += fmtVal(lst[i].max,n,d)
    print(msg)
    print(hdrmods1)
    print(hdrmods2)
    print(ln1)
    print(ln2)
    print(ln3)
    print(ln4)


printPMTStats('NPE Statistics:',s_npe_pmt)
printPMTStats('Ratio Statistics:',s_rat_pmt,d=2)


## Finally, print the module ratios
print()
print('Module Normalizations by PMT')
print(hdrmods1 + '    AVG       STDV       MIN       MAX')
print(hdrmods2 + '   ------    ------    ------    ------')
s_rat2_mod = []
for im,m in enumerate(modules):
    s_rat2_mod.append(StatList())
    ln1 = '[%02d-%02d]'%(m.row,m.col)
    for ipmt,pmt in enumerate(m.pmts):
        rv = s_npe_pmt[ipmt].rvals[im]
        v  = s_npe_pmt[ipmt].vals[im]
        ln1 += fmtVal(rv,d=3,hi=(v==s_npe_pmt[ipmt].max))
        s_rat2_mod[im].append(s_npe_pmt[ipmt].rvals[im],1.0)
    s_rat2_mod[im].Finalize()
    ln1 += fmtVal(s_rat2_mod[im].avg,d=3)
    ln1 += fmtVal(s_rat2_mod[im].std,d=3)
    ln1 += fmtVal(s_rat2_mod[im].min,d=3)
    ln1 += fmtVal(s_rat2_mod[im].max,d=3)
    print(ln1)
