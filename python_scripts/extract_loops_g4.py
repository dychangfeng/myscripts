import pandas as pd
import re
import numpy as np

hg38_all = pd.read_table('/uufs/chpc.utah.edu/common/home/u0674686/ncbi_genomes/bed_files_for_all_genomes/all_G4/hg38_all.fa.all_G4.bed', header = None)


def reverse_comp(seq):
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join(bases_dict.get(base,base) for base in seq[::-1])

print tt, reverse_comp(tt)

l1 = []
l2 = []
l3 = []
l4 = []
l5 = []
gx_runs_g4 = []
ls = []
gs = re.compile(r'([gG]{3,})')
for j in range(len(hg38_all.iloc[:,6])):
    tt = hg38_all.iloc[:,6].iloc[j]
    starts = []
    ends = []
    if tt[0] == 'G' or tt[0] == 'g': 
        tt = tt
    else :
        tt = reverse_comp(tt) ## convert to G4s instead of C4s
    s = ''
    for m in gs.finditer(tt):
        starts.append(m.start())
        ends.append(m.end())
        s = s + str(len(m.group())) + ':'
    if len(starts)<4:  ## ingore the G4s with G loops only, or ignore the polyG sequence
        gx_runs_g4.append(tt)
        continue
    ls.append(s)
    if len(starts)==5:
        l5.append('None')
    elif len(starts)==4:
        l4.append('None')  
        l5.append('None')
    for i in range(len(starts)-1):
        if i == 0:
            l1.append(tt[ends[i]:starts[i+1]]) ## the first start is not a loop, it is a G run, so it is i + 1
        elif i == 1:
            l2.append(tt[ends[i]:starts[i+1]])
        elif i == 2:
            l3.append(tt[ends[i]:starts[i+1]])
        elif i == 3:
            l4.append(tt[ends[i]:starts[i+1]])
        elif i == 4:
            l5.append(tt[ends[i]:starts[i+1]])
hg38_loops_pd = pd.DataFrame({'loop1':l1,'loop2':l2,'loop3':l3,'loop4':l4,'loop5':l5, 'g_runs':ls})
    

print len(l1), len(l2), len(l3), len(l4), len(l5), len(gx_runs_g4)



hg38_loops_pd.to_csv('/uufs/chpc.utah.edu/common/home/u0674686/ncbi_genomes/hg38_all_g4_loops.csv')

