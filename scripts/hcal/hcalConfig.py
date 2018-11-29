#!/usr/bin/env python3

import re

maxmaxnpe=7806.0 ## This is the PMT with the highest QE for reference


def lineCleanup(ln):
    ln = ln.rstrip() ## Remove trailing new line characters
    ln = ln.split('#')[0] ## Remove any comments
    return ln

class SummaryValue(object):
    __slots__ = ('row','col','led','npe','err','sig','n','idx')
    def __init__(self,items=[]):
        if len(items) == 7:
            self.row = int(items[0])
            self.col = int(items[1])
            self.led = int(items[2])
            self.npe = float(items[3])
            self.err = float(items[4])
            self.sig = float(items[5])
            self.n   = int(items[6])
            self.idx = self.getidx(self.row,self.col,self.led)
    def __lt__(self,other):
        return self.idx < other.idx

    @staticmethod
    def getidx(row,col,led):
        return led*1000+(row-1)*12+(col-1)

def readSummaryFile(run):
    f = open('summary/summary_%d.dat'%(run),'r')
    vals=[]
    for l in f:
        l = l.rstrip();
        items=l.split();
        if len(items) == 7:
            vals.append(SummaryValue(items))
            #if int(items[2]) == led and int(items[1]) == col:
            #    vals[int(items[0])] = float(items[3])
    ## Sort the values on return (for easier processing later)
    return sorted(vals)

def getSummary(summary,row,col):
    return next((sv for sv in summary if sv.row == row and sv.col == col and sv.led==32),None)


def readHCalList(filename='configs/pmt.cfg',by_run=True,read_summary=False,unique=True):
    current_run = -1
    f = open(filename,'r')
    lc = 0
    DB={}
    summary = []
    for l in f:
        lc += 1
        l = lineCleanup(l)
        if not l: ## Skip empty lines
            continue
        ## Check if a new run is specified
        m = re.match(r'\[(\d*)\]',l)
        if m:
            current_run=int(m.group(1))
            if by_run:
                DB[current_run] = []
            if current_run>0 and read_summary:
                summary = list(readSummaryFile(current_run))
        else:## Should be list of items
            items = l.split()
            litems = len(items)
            if litems < 3:
                continue
            row=int(items[0])
            col=int(items[1])
            pmt=items[2]
            extra = None
            if len(items) > 3:
                extra=items[3]
            sv = None
            if read_summary:
                sv = getSummary(summary,row,col)
                if sv == None:
                    print('Run %4d [%2d-%2d] %7s'%(current_run,row,col,pmt))
            if by_run:
                DB[current_run].append([int(row),int(col),pmt,sv])
            elif extra == None:
                if unique and pmt in DB:
                    print('Line %5d : PMT %s in run %d already in run %d'%(lc,pmt,current_run,DB[pmt][2]))
                else:
                    DB[pmt] = [int(row),int(col),current_run,sv]

            #if store_by_pmt and unique and key in DB:
            #else:
            #    DB[key] = [int(row),int(col),val]
    return DB

def readRunList(filename='configs/pmt.cfg',read_summary=False,unique=False):
    return readHCalList(filename,True,read_summary,unique)

def readPMTList(filename='configs/pmt.cfg',read_summary=True,unique=True):
    unsorted_pmts=readHCalList(filename,False,read_summary,unique)
    sorted_pmts={}
    for pmt in sorted(unsorted_pmts.keys()):
        sorted_pmts[pmt] = unsorted_pmts[pmt]
    #return list(sorted(readHCalList(filename,False,read_summary).keys()))
    ## Sort them for faster access later
    #return pmts
    return sorted_pmts


def readNorms(cfg='configs/normvals.cfg'):
    f = open(cfg,'r')
    norms=[ [-1.for c in range(12)] for r in range(24) ]
    for l in f:
        l = lineCleanup(l)
        if not l: ## Skip empty lines
            continue
        items = l.split()
        if len(items) != 3:
            continue
        row=int(items[0])
        col=int(items[1])
        norm=float(items[2])
        if row > 0 and row < 24 and col > 0 and row < 12:
            norms[row][col] = norm

    return norms

