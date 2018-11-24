# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 14:27:25 2017

@author: Yun
"""
from collections import Counter
import pandas as pd
import re
from IPython.display import HTML, display

def reverse_comp(seq):
    """take a sequence and return the reverse complementary strand of this seq
    input: a DNA sequence
    output: the reverse complementary of this seq"""
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join(bases_dict.get(base,base) for base in seq[::-1])

def get_G4_list_from_df_bed(hg38_all):
    """take a dataframe from bedfile of Quadpaser, convert C4 to G4s, return all G4 sequences as a list"""
    g4_ls =[]
    #gs = re.compile(r'([gG]{3,})')
    for j in range(len(hg38_all.iloc[:,6])):
        tt = hg38_all.iloc[:,6].iloc[j]
        tt =tt.upper() ## convert to upper case
        if tt[0] == 'G':
            tt = tt
        else :
            tt = reverse_comp(tt) ## convert to G4s instead of C4s
        g4_ls.append(tt)
    return g4_ls

def count_most_frequent_G4_from_G4list(g4ls, n=15):
    """find the most frequent G4 from a list of g4s,
    return a dataframe with 3 columns['G4', 'Count', 'G4_length']"""
    g4_counter = Counter(g4ls) # convert list to Counter object
    g4_most = g4_counter.most_common(n) ## figure out the most common 15 sequences
    df = pd.DataFrame.from_records(g4_most, columns=['G4', 'Count']) ## convert to dataframe, rename column name
    df['G4_length']=[len(x) for x in df.G4] ## add the length of G4 as the third column
    return df

def highlight_G_runs(df):

    """highlight the G runs in the top G4s to compare the loops
    take dataframe from 'count_most_frequent_G4_from_G4list'
    display highlighted G runs in G4s"""

    regex=re.compile("[G]{3,}")
    i = 0
    for text in df.G4:
        output=''
        for m in regex.finditer(text):
            output += "".join([text[i:m.start()],'<font color="red">%s</font>'%text[m.start():m.end()]])
            i = m.end()
        display(HTML("".join([output, text[m.end():]])))
