# -*- coding: utf-8 -*-
import pandas as pd
import re
import numpy as np
def reverse_comp(seq):
    """take a sequence and return the reverse complementary strand of this seq
    input: a DNA sequence
    output: the reverse complementary of this seq"""
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join(bases_dict.get(base,base) for base in seq[::-1])

def get_5_loops_from_bedfile(df):
    """input: dataframe of the bedfile from quadrpaser
    Only consider the first 5 loops, for G4s have less than 5 loop, fill with np.nan
    output: a new dataframe for all the loops for each G4"""

    l1 = []
    l2 = []
    l3 = []
    l4 = []
    l5 = []
    #gx_runs_g4 = []
    ls = []
    gs = re.compile(r'([gG]{3,})')
    for j in range(len(df.iloc[:,6])):
        tt = df.iloc[:,6].iloc[j].upper()
        starts = []
        ends = []
        if tt[0] == 'G':
            tt = tt
        else :
            tt = reverse_comp(tt) ## convert to G4s instead of C4s
        s = ''
        for m in gs.finditer(tt):
            starts.append(m.start())
            ends.append(m.end())
            s = s + str(len(m.group())) + ':'
        if len(starts)<4:  ## ingore the G4s with G loops only, or ignore the polyG sequence
            #gx_runs_g4.append(tt)
            continue
        ls.append(s) ## the type of G runs, eg, 3:3:3:3:
        if len(starts)==5:
            l5.append(np.nan)
        elif len(starts)==4:
            l4.append(np.nan)
            l5.append(np.nan)
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
    return pd.DataFrame({'loop1':l1,'loop2':l2,'loop3':l3,'loop4':l4,'loop5':l5, 'g_runs':ls})
