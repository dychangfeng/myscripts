import os
import sys
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd

"""get .fa file for the sequences flanking the TSS, read in one by one and calculated the GC for defined bin size such as 200bp """ 
fa_file ='/Users/Yun/Documents/human_cruciform/sequence_around_TSS/sequence_aroundTSS_4000_strand.fa'
def get_seq_dict(fa_file):
    seq_dict = {} # dic to save all the GC for each sequence
    with open(fa_file, 'r') as fh:
        lines=fh.readlines()
        name=''
        for line in lines:
            if line.startswith('>'):
                name=line
                seq_dict[name]='' #otherwise key error
            else:
                seq_dict[name] +=line
    return(seq_dict)
def get_GC(seq, base='GC'):
    """caculate GC percent in a sequence
    base can be 'G', 'C' or 'CG'
    by default, base='GC', which calculate the average GC content"""
    seq=seq.upper()
    if base=='GC':
         GC=(seq.count('G')+seq.count('C'))*100/len(seq) #
    elif base == 'G':
         GC=(seq.count('G'))*100/len(seq) 
    elif base == 'C':
         GC=(seq.count('C'))*100/len(seq)  
    return(int(round(GC)))

def GC_in_bins(seq, bin_width = 100, base='GC'):
    seqGC=np.ones(int(4000/bin_width), dtype=int)
    for i in range(0,4000,bin_width):
         seqGC[int(i/bin_width)] = get_GC(seq[i:i+bin_width], base=base)
    return(seqGC)

def GC_dict(sequence_dict, bins=40, base='GC'):
    """calcuate GC for a dictionary of sequences"""
    n=len(sequence_dict)
    GCcounter=np.ndarray((n, bins), dtype=int)
    My_keys= sorted(sequence_dict.keys())
    for i in range(n):
        GCcounter[i,:] = GC_in_bins(sequence_dict[My_keys[i]], bin_width=int(4000/bins), base=base) ## use map function there??
    return(GCcounter)
bins = 200

myArr=GC_dict(get_seq_dict(fa_file),bins=bins)
print(myArr.mean(axis=0))
#tt=pd.DataFrame(myArr, columns=range(-19,20,2))
allMeans = myArr.mean(axis=0)
plt.scatter(range(-1990,2000, 20), allMeans)
#plt.ylim(10,40)
plt.savefig('/Users/Yun/Documents/human_cruciform/sequence_around_TSS/means_C20.png')
# ==============================================================================
# calucate average GC, G, or C for this file
# ==============================================================================
all_mean = np.ndarray((bins,4), dtype=float)
all_mean[:,0] = range(-1999,2000, 20)
bases=['G', 'C', 'GC']
for i in range(3):
     myArr=GC_dict(get_seq_dict(fa_file),bins=bins, base=bases[i]) 
     all_mean[:,i+1] = myArr.mean(axis=0)
col_head = ['bins'] + bases
df = pd.DataFrame(all_mean, columns=col_head)
print(df.head())
df.to_csv('/Users/Yun/Documents/human_cruciform/sequence_around_TSS/mean_GC_G_20size.csv', sep=',') ## save to file as float, recover using function fromfile()
