
# load packages

import os,sys
from numpy import *
#import PyQt5 # Only when not anaconda2 is not sourced
#from matplotlib import *
#import scipy as scp

# Load File function
def load_FILE(fname,spt='\t'):

    infile = open(fname,'r')
    lines = infile.readlines()
    peaks = []
    for line in lines:
        cols = line.split(spt)
        cols[-1] = cols[-1].strip()
        peaks.append(cols)
    peaks = array(peaks)

    return peaks

## Load diff-exp table from R
tbl = load_FILE('Voom_diffexp_gnms.txt') # cutoff=1 after normalization
hdr = tbl[0]
tbl = tbl[1:]
gnms = tbl[:,0]
ensid = tbl[:,1]
fc = tbl[:,2].astype(double)
avgx = tbl[:,3].astype(double)
pval = tbl[:,5].astype(double)
qval = tbl[:,6].astype(double)
mlogp = -log10(pval)
mlogq = -log10(qval)

###########################################################################

## up and downregulated genes

# Alternative 1: cutoff based on q-value

ctf = 0.05
id0 = where(qval <= 0.05)[0]
id01 = where(fc > 0)[0]
id02 = where(fc < 0)[0]
id1 = intersect1d(id0,id01)
id2 = intersect1d(id0,id02)
str1 = 'q > '+str(ctf)
str2 = 'q < '+str(ctf)
idd = union1d(id1,id2)

# Alternative 2: Cutoff based on top NN genes
NN = 100
id0 = arange(NN)
id01 = where(fc > 0)[0]
id02 = where(fc < 0)[0]
id1 = intersect1d(id0,id01)
id2 = intersect1d(id0,id02)
str1 = 'top '+str(NN)+' up'
str2 = 'top '+str(NN)+' down'
idd = union1d(id1,id2)

############################################################################

## Lecture plot
str1 = 'Up'
str2 = 'Down'

## MA plot

ff = figure()
ff.set_figwidth(12)
ff.set_figheight(10)

plot(avgx,fc,'k.',ms=2)
plot(avgx[id1],fc[id1],'r.',ms=10,label=str1)
plot(avgx[id2],fc[id2],'b.',ms=10,label=str2)

xlabel('Average log-expression',size=20)
ylabel('log-fold-change',size=20)

fsz = mpl.font_manager.FontProperties(size=20) 
legend(prop=fsz)

plt.savefig('MA_plot.pdf',format='pdf')

#############################################################################

## Volcano Plot - p-value

ff = figure()
ff.set_figwidth(12)
ff.set_figheight(10)

plot(fc,mlogp,'k.',ms=2)
plot(fc[id1],mlogp[id1],'r.',ms=10,label=str1)
plot(fc[id2],mlogp[id2],'b.',ms=10,label=str2)

xlabel('Log2 Fold Change',size=20)
ylabel('-log10(P-Value)',size=20)

fsz = matplotlib.font_manager.FontProperties(size=20) 
legend(prop=fsz)

## add gene-names
for i in idd:
    tt = text(fc[i],mlogp[i],gnms[i],fontsize=8)

plt.savefig('Volcano_plot.pdf',format='pdf')

#############################################################################
#############################################################################
#############################################################################
