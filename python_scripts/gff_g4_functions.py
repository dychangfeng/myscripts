"""modules Description
Copyright (c) 2017 Yun Ding <u0674686@utah.edu>
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@author:  Yun Ding
@contact: u0674686@utah.edu
This module contains functions to calculate G4s in cds, around TSS regions, and UTR regions"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pybedtools

def reverse_comp(seq):
    """take a sequence and return the reverse complementary strand of this seq
    input: a DNA sequence
    output: the reverse complementary of this seq"""
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    ##get the complementary if key is in the dict, otherwise return the base
    return ''.join(bases_dict.get(base, base) for base in seq[::-1])

def g4_overlap_cds(base_dir):
    """get base direction of all_G4_clean and gff_files,
    calculate the percentage of G4 overlap with CDS, and CDS percentage
    return a data frame with three columns['GCA_access', 'G4_in_CDS', 'cds_Mb']"""
    #base_dir='/Users/Yun/Documents/bacteria_G4/D_thermus/'
    base_dir_gff = base_dir + 'gff_files/'
    cds_intersect_g4 = []
    for i in os.listdir(base_dir_gff):
        gff = pybedtools.BedTool(os.path.join(base_dir_gff+i))
        try:
            g4_bed = pybedtools.BedTool(\
            os.path.join(base_dir+'all_G4_clean/'+i.split('.gff.gz')[0]+'.fna.all_G4.bed'))
        except:
            continue
        cds_gff = gff.filter(lambda x: str(x[2]) == 'CDS')
        ## after filter you can only use it once because it is generator
        cds_percentage = len(cds_gff.intersect(g4_bed))/float(len(g4_bed))
        ## fraction of g4 intersect with cds
        cds_gff = gff.filter(lambda x: str(x[2]) == 'CDS').saveas(base_dir+'_cds.txt').sort().merge()
        ## regenerate the generator for the following for loop
        ## need to consider overlap problem later??
        s_cds = 0 ## use s to calculate the length of CDS
        for cds_range in cds_gff:
            s_cds += len(cds_range)
        cds_intersect_g4.append(('_'.join([i.split('_')[0], i.split('_')[1]]), \
        cds_percentage, float(s_cds)/1000000))
        ## cds in Million basepaires
    cds_g4 = pd.DataFrame.from_records(cds_intersect_g4)
    cds_g4.rename(columns={0:'GCA_access', 1: 'G4_in_CDS', 2:'cds_Mb'}, inplace=True)
    return cds_g4

def g4_tss_window(base_dir, window=100):
    """calculate how many g4 within each window of around TSS
    the defaut window is 100
    save the file to disk"""
    #base_dir='/Users/Yun/Documents/bacteria_G4/D_thermus/'
    tss_dir = base_dir + 'tss_files/'
    distr_g4 = {}
    ## create a dict with the file name as the key and the list of number of G4s
    # within the window as the value
    for i in os.listdir(tss_dir):
        distr = []
        tss = pybedtools.BedTool(os.path.join(tss_dir+i))
        try:
            g4_bed = pybedtools.BedTool(os.path.join(\
            base_dir+'all_G4_clean/'+''.join(\
            [i.split('.bed')[0], '.fna.all_G4.bed'])))
        except:
            continue # the bed file might not exist
        for window_tss in range(0, window, 10):
            g4_in_range = g4_bed.window(tss, l=window_tss, r=window_tss, u=True)
            n_g4 = float(len(g4_bed))
            distr.append(len(g4_in_range)/n_g4)
        distr_g4[i.split('.bed')[0]] = distr
    df1 = pd.DataFrame.from_records(distr_g4)
    df1['row_name'] = pd.Series(range(0, window, 10))
    df1.set_index('row_name', inplace=True)
    df2 = df1.transpose()
    df2.reset_index(inplace=True)
    df2['index'] = pd.Series(\
    ['_'.join([x.split('_')[0], x.split('_')[1]]) for x in df2['index']])
                                  ## clean up the column with all the file name
    df2.rename(columns={'index':'GCA_access'}, inplace=True)
    df2.to_csv(''.join([base_dir, "G4_within_window_of_TSS.txt"]), \
    sep='\t', index=False)
    ## save to a new file

def plot_g4_tss2(ls_distance, start=-300, end=300, step=30, title='', base_dir=''):
    """plot, histgram for a list of numbers
    input: list of distance between G4 and TSS
    save png figure in the figures folder in the base directory"""
    #axes = plt.gca() # get current axes
    fig, axes = plt.subplots()
    font = {'family':'serif', 'color':'darkred', 'weight':'normal', 'size':20}
    axes.set_xlim([start, end])
    axes.set_xlabel('Distance (nt)', fontdict=font)
    axes.set_ylabel('Counts', fontdict=font)
    axes.set_title(title, fontdict=font)
    bins = np.linspace(start, end, (end-start)/step)
    plt.hist(ls_distance, bins=bins, color='red', alpha=0.5)
    matplotlib.rc('xtick', labelsize=16) ## run after the plot function
    matplotlib.rc('ytick', labelsize=16)
    plt.tight_layout()
    fig.savefig(os.path.join(base_dir + 'figures/' + title + '.png'), dpi = 200)


def get_g4_distance_to_closest_tss_3(base_dir, GCA_list):
    """in the base directory, , if the GCA is in GCA _list, calculate the the distance of G4 to the closest
returns  dict with GCA number as the key as distance list as teh value"""
    counts_around_tss = {}
    for i in os.listdir(base_dir+'tss_files/'):
        name = i.split('_')
        gca = '_'.join(name[0:2])
        if gca in GCA_list:
            tss = pybedtools.BedTool(base_dir+'tss_files/'+i)
            i = i.split('.bed')
            ## if there are G4 in the genes in the corresponding genome
            if os.path.isfile(base_dir+'all_G4_clean/' + i[0] + '.fna.all_G4.bed'):
                g4_bed = pybedtools.BedTool(base_dir+'all_G4_clean/' + i[0] + \
                '.fna.all_G4.bed')
                if len(g4_bed) > 0:
                    ## 'b' means tss strandness is used to determin upstream and downstream
                    g4_closest_tss = pd.DataFrame.from_records(g4_bed.closest(tss, D='b'))
                    g4_closest_tss=g4_closest_tss.ix[g4_closest_tss[6]!=u'.',:]
                else:
                    pass
        ## get the GCA number of a file as the key, the distance is the values
            counts_around_tss[gca] = [int(t) for t in g4_closest_tss.iloc[:, -1]]
        else:
            pass

    return counts_around_tss

def get_g4_distance_to_closest_tss(base_dir, G4_dir = 'all_G4_clean/', pos_dir = 'tss_files/', start=-100, end=100, bins=20, strand = "both"):
    """in the base directory, calculate the the distance of G4 to the closest TSS in each file
    pos_dir is the directory to look for tss. 
    It can be 'gene_end_files/' when examing enrichment in 3UTR regions
    start and end are the distance around TSS, bins is the number of bins
    return a data frame with all the distance and fold enrichment"""
    counts_around_tss = {}
    G4_tail = '_'.join(G4_dir.split('_')[0:2])
    for i in os.listdir(base_dir+pos_dir):
        tss = pybedtools.BedTool(base_dir+pos_dir+i)
        i = i.split('.bed')
        ## if there are G4 in the genes in the corresponding genome
        G4_path = base_dir+G4_dir+ i[0] + '.fna.' + G4_tail+ '.bed' 
        if os.path.isfile(G4_path):
            g4_bed = pybedtools.BedTool(G4_path)
            if len(g4_bed) > 0:
                ## 'b' means tss strandness is used to determin upstream and downstream
                g4_closest_tss = pd.DataFrame.from_records(g4_bed.closest(tss, D='b'))
                ## remove rows without annotation in the tss file, the gene name is '.',
                #and all other values are '-1'
                g4_closest_tss=g4_closest_tss.ix[g4_closest_tss[6]!=u'.',:]
                if strand == "template": # template strand only
                	g4_closest_tss = g4_closest_tss.ix[g4_closest_tss[5]!=g4_closest_tss[11],:]
                elif strand == "code":
                	g4_closest_tss = g4_closest_tss.ix[g4_closest_tss[5]==g4_closest_tss[11],:]
            else:
                pass
            ## get the last column (distance) into bins and counts
            counts, bin_edge = np.histogram(\
            [int(t) for t in g4_closest_tss.iloc[:, -1]], \
            bins=np.linspace(start, end, bins+1))
            name = i[0].split('_')
            counts_around_tss[name[0] + '_' + name[1]] = counts
        else:
            continue
    g4_counts_around_tss = pd.DataFrame.from_records(counts_around_tss).transpose()
    g4_counts_around_tss['G4_around_tss_mean'] = g4_counts_around_tss.apply(\
    np.mean, axis=1)
    ## add the mean counts of G4 around TSS
    g4_counts_around_tss['fold_enrichment'] = g4_counts_around_tss[bins/2]/g4_counts_around_tss['G4_around_tss_mean']
    ## calculate the fold enrichment of G4_tss over the mean
    g4_counts_around_tss.reset_index(inplace=True)
    # the center bin is the TSS which is : bins/2 -1
    g4_counts_around_tss.rename(columns=\
    {'index':'GCA_access', bins/2:'G4_counts_at_tss'}, inplace=True)
    return g4_counts_around_tss

def get_g4_distance_to_closest_tss_2(base_dir, strand = "both"):
    """in the base directory, calculate the the distance of G4 to the closest
returns  dict with GCA number as the key as distance list as teh value,
strand is an option to choose either code for the coding, or template for the template strand"""
    counts_around_tss = {}
    for i in os.listdir(base_dir+'tss_files/'):
        tss = pybedtools.BedTool(base_dir+'tss_files/'+i)
        i = i.split('.bed')
        ## if there are G4 in the genes in the corresponding genome
        if os.path.isfile(base_dir+'all_G4_clean/' + i[0] + '.fna.all_G4.bed'):
            g4_bed = pybedtools.BedTool(base_dir+'all_G4_clean/' + i[0] + \
            '.fna.all_G4.bed')
            if len(g4_bed) > 0:
                ## 'b' means tss strandness is used to determin upstream and downstream
                g4_closest_tss = pd.DataFrame.from_records(g4_bed.closest(tss, D='b'))
                g4_closest_tss=g4_closest_tss.ix[g4_closest_tss[6]!=u'.',:]
                if strand == "template": # template strand only
                	g4_closest_tss = g4_closest_tss.ix[g4_closest_tss[5]!=g4_closest_tss[11],:]
                elif strand == "code":
                	g4_closest_tss = g4_closest_tss.ix[g4_closest_tss[5]==g4_closest_tss[11],:]
            else:
                pass
        name = i[0].split('_')
        ## get the GCA number of a file as the key, the distance is the values
        counts_around_tss['_'.join([name[0],  name[1]])] = [int(t) for t in g4_closest_tss.iloc[:, -1]]

    return counts_around_tss

def get_end_position_from_gff(file_name, base_dir):
    """take the file name and base_dir, return the file with end position of
    as bed file format"""
    ucsc_tss=[]
    with open(base_dir+file_name, 'r') as f0:
        lines=f0.readlines()
        for line in lines:
            line=line.split('\t')
            if line[6]=='-':
                ucsc_tss.append((line[0], line[3], line[3], line[5], line[5], line[6]))
            elif line[6]=='+':
                ucsc_tss.append((line[0], line[4], line[4], line[5], line[5], line[6]))
    with open(base_dir+file_name+'.end.bed', 'w') as f0:
        f0.write('\n'.join('{}\t{}\t{}\t{}\t{}\t{}'.format(\
                        x[0], x[1], x[2], x[3], x[4], x[5]) for x in ucsc_tss))

def most_common_G4_around_tss(base_dir, G4_dir = 'all_G4_with_seq/', n = 30):
    """in the base directory, calculate the the distance of G4 to the closest
    TSS in each file
    find the most common G4 that within 30 bps of tss
    n: the number of bases away from TSS"""
    G4_around_tss = {}
    for i in os.listdir(base_dir+'tss_files/'):
        tss = pybedtools.BedTool(base_dir+'tss_files/'+i)
        i = i.split('.bed')
        ## if there are G4 in the genes in the corresponding genome
        if os.path.isfile(base_dir+G4_dir+ i[0] + '.fna.' + ''.join(G4_dir.split('_')[0:2]) + '.bed'):
            g4_bed = pybedtools.BedTool(base_dir+ G4_dir + i[0] + \
            '.fna.all_G4.bed')
            if len(g4_bed) > 0:
                ## 'b' means tss strandness is used to determin upstream and downstream
                g4_closest_tss = pd.DataFrame.from_records(g4_bed.closest(tss, D='b'))
                ## remove rows without annotation in the tss file, the gene name is '.',
                #and all other values are '-1'
                g4_closest_tss=g4_closest_tss.ix[g4_closest_tss[6]!=u'.',:]
                g4_closest_tss['dist']=[int(j) for j in g4_closest_tss.iloc[:,-1]]
                g4_tss = g4_closest_tss.loc[abs(g4_closest_tss.dist)<n,]
                g4_ls = []
                for t in range(g4_tss.shape[0]):
                    seq = g4_tss.iloc[t,3]
                    seq = seq.encode('ascii',  'ignore')
                    if seq[0].upper() == 'G':
                        g4_ls.append(seq)
                    else:
                        g4_ls.append(reverse_comp(seq))

            else:## pass when there is no G4
                pass
        name = i[0].split('_')
        G4_around_tss[name[0] + '_' + name[1]]=Counter(g4_ls).most_common(1) ## count the list
    return G4_around_tss
