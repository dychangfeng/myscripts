# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 14:27:25 2017

@author: Yun
"""
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

def reverse_comp(seq):
    """take a sequence and return the reverse complementary strand of this seq
    input: a DNA sequence
    output: the reverse complementary of this seq"""
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join(bases_dict.get(base,base) for base in seq[::-1])

def get_single_type_G4(df, g4_type='3:3:3:3:'):
    """input: a dataframe from extract_loops.get_5_loops_from_bedfile
    g4_type with defaut 3:3:3:3:
    return one type of G4 as a dataframe with sorted loop frequence and length
    return df: value(loop counts), loop, loop_length"""
    if g4_type!='None':
        df= (df[df.g_runs == g4_type]).iloc[:,0:int(len(g4_type)/2)]
    else:
        df=df ## return melted all loops when g4_type="None"
    melt_df= pd.melt(df, id_vars = ['g_runs'], var_name = 'all_loop')
    all_loop_df=(melt_df.value.value_counts()).to_frame()
    all_loop_df['loop']=all_loop_df.index
    all_loop_df['loop_length'] = [len(x) for x in all_loop_df.loop]
    return all_loop_df

def plot_loop_for_single_typeG4(df, n=20):
    """input: a melted dataframe from get_single_type_G4
    output: a barplot with top n loops """
    df.sort_values(by='value', ascending=False, inplace=True)
    sns.set_style('dark') ## set themes
    sns.set_style('ticks')
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 8)
    plt.xticks(rotation=60)
    sns.set_style('ticks')
    sns.set_context('talk')
    sns.set(font_scale=1)
    colors = sns.color_palette("husl", n_colors=len(set(df.loop_length[0:n])))
    color_dict = dict(zip(set(df.loop_length[0:n]), colors))
    palette = (df[0:n])['loop_length'].map(color_dict)
    g = sns.barplot(data = df[0:n], x = 'loop', y = 'value', palette = palette)
    g.set_xlabel("Loops", fontsize=20)
    g.set_ylabel("Counts of each specific loop", fontsize = 20)
    g.tick_params(labelsize=16)
    plt.tight_layout()
    ##---add label-----------------------
    handles=[]
    for i in range(len(colors)):
        handles.append(mpatches.Patch(color=colors[i], label=' '.join([str(i+1), 'nt loop'])))
    plt.legend(handles=handles, prop={'size':16})
