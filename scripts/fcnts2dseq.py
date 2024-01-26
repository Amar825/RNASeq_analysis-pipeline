# converts gene counts from featureCounts to table format usable for DESeq
import os,sys
from numpy import *

def write_FILE(table,fname,delim='\t'):

        outfile = open(fname,'w')
        for line in table:
                outfile.write(delim.join([str(s) for s in line])+os.linesep)
        outfile.close()

def load_FILE(fname,spt='\t'):

    infile = open(fname,'r')
    lines = infile.readlines()
    lines = lines[1:]
    peaks = []
    for line in lines:
        cols = line.split(spt)
        cols[-1] = cols[-1].strip()
        ncols = hstack((cols[0],cols[6:]))
        peaks.append(ncols)
    peaks = array(peaks)

    return peaks

fname = sys.argv[1]
peaks = load_FILE(fname)
fname_out = fname[:-4]+'.tsv'
write_FILE(peaks,fname_out)
