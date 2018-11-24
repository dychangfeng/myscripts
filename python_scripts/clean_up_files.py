
"""modules Description
Copyright (c) 2017 Yun Ding <u0674686@utah.edu>
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@author:  Yun Ding
@contact: u0674686@utah.edu
This module clean up the G4_bed files, find tss from GFF3 files,
and and the overlap of G4 with coding regions and TSS regions
fasta_parser is a function to parse fasta files into dict"""
import os
import gzip

def get_start_position_from_gff(file_name, base_dir):
    """take the file name and base_dir which is the gff_files folder, return the file with start position of
    as bed file format
    line[8].split()[0] is the ID"""
    ucsc_tss=[]
    with open(base_dir+file_name, 'r') as f0:
        lines=f0.readlines()
        for line in lines:
            line=line.split('\t')
            if len(line[0])>5: ## ignore sequences not in chromosome
                continue
            if line[0].startswith('#'):
                continue
            elif line[6]=='+':
                ucsc_tss.append((line[0], line[3], line[3], line[5], line[8].split(';')[0], line[6]))
            elif line[6]=='-':
                ucsc_tss.append((line[0], line[4], line[4], line[5], line[8].split(';')[0], line[6]))
    with open(base_dir+file_name+'.bed', 'w') as f0:
        f0.write('\n'.join('{}\t{}\t{}\t{}\t{}\t{}'.format(\
                        x[0], x[1], x[2], x[3], x[4], x[5]) for x in ucsc_tss))

def get_tss_list(file_path, pos1 = 3, pos2 = 4):
    """take a gff3 file object
    return a list of tuples with all the genes from gff3 file
        line[8].split()[0] is the ID, add gene ID so later can extract the gene ID with distance
    To Do: add lines so it can open both .gz files and flat files
    for TSS start: pos1 = 3, pos2 = 4
    for the end of transcription site: pos1 = 4, pos2 = 3
    The default is looking for TSS
    line[8].split(';')[1] is the gene symbol
"""
    if '.gz' in file_path:
        with gzip.open(file_path) as f0:
            lines = f0.readlines()
            genes = []
            for line in lines:
                if not line.startswith('#'):
                    line = line.split('\t')
                    if line[2] == 'gene': ## only pick out 'genes'from all the features
                        if line[6] == '+':
                            genes.append((line[0], int(line[pos1]), \
                            int(line[pos1]), '.', line[8].split(';')[1].rstrip(), line[6]))
                            ## append seq_id, tss_start(twice), strand
                        else:
                            genes.append((line[0], int(line[pos2]), \
                            int(line[pos2]), '.', line[8].split(';')[1].rstrip(), line[6]))
                        # on '-' strand, append the end, instead of start
    if not '.gz' in file_path: # not a gzip files
        with open(file_path) as f0:
            lines = f0.readlines()
            genes = []
            for line in lines:
                line = line.split('\t')
                if not line[0].startswith('#'):
                    if line[2] == 'gene': ## only pick out 'genes'from all the features
                        if line[6] == '+':
                            genes.append((line[0], int(line[pos1]), \
                            int(line[pos1]), '.', line[8].split(';')[1].rstrip(), line[6]))
                            ## append seq_id, tss_start(twice), strand
                        else:
                            genes.append((line[0], int(line[pos2]), \
                            int(line[pos2]), '.', line[8].split(';')[1].rstrip(), line[6]))

    if len(genes) > 0:
        return sorted(genes, key=lambda x: x[0:2])
        ## sort by gene_seq, start, and end

    else: 
        return "no genes in this file"

def create_tss_from_gff3(base_dir, start = True):
    """input: base directionary
    output: save tss as bed6 format
    to do: make sure not to include to .DS_store file during analysis
    the default is looking for the tss
    To look for the end of a gene: start = False, pos1 = 4, pos2 =3 """
    #base_dir='/Users/Yun/Documents/bacteria_G4/D_thermus/'
    #base_dir_gff = base_dir + 'gff_files/'
    if start==True: #check whether get the beginning or the end of the file
        store_dir =base_dir + 'tss_files/'
        pos1 = 3
        pos2 = 4
        if not os.path.isdir(base_dir + 'tss_files/'):
            os.mkdir(base_dir + 'tss_files/')
    else: ## create a folder for the end of the gene sites
        store_dir =base_dir + 'gene_end_files/' 
        pos1 = 4
        pos2 = 3
        if not os.path.isdir(base_dir + 'gene_end_files/'):
            os.mkdir(base_dir + 'gene_end_files/') 
                       
    for i in os.listdir(base_dir + 'gff_files/'):
        if i.startswith('.'):
            continue ## ingore hidden files such as .DS_store
        elif os.path.isfile(store_dir + i):
            print("file {} exist".format(i))
            continue
        else:
            gff = get_tss_list(base_dir + 'gff_files/'+i, pos1=pos1, pos2=pos2)
            ## get only genes out of gff3 files
            if not type(gff) == str:
                ## when there is no genes in gff3, read_gff3 return a string
                with open(os.path.join(store_dir, \
                i.split('.gff')[0] + '.bed'), 'w') as f0:
                    f0.write('\n'.join(['{}\t{}\t{}\t{}\t{}\t{}'.format(\
                    x[0], x[1], x[2], x[3], x[4], x[5]) for x in gff]))
            else:
                print("no genes in file: %s"%(i))

def trim_g4_chr(base_dir, G4_dir = 'all_G4/', G4_clean_dir = 'all_G4_clean/'):
    """trim off the extra text to match chromsome name in
    the GFF files to use bedfiles,
    save the new file to the folder all_G4_clean"""
    #base_dir='/Users/Yun/Documents/bacteria_G4/D_thermus/'
    G4_dir = base_dir + G4_dir
    if not os.path.isdir(base_dir + G4_clean_dir):
        os.mkdir(base_dir + G4_clean_dir)
    for i in os.listdir(G4_dir):
        if i.startswith('.'):
            continue # igore hiden files such as .DS_store
        with open(G4_dir+i, 'r') as fp:
            lines = fp.readlines()
            newlines = []
            for line in lines:
                line = line.split('\t')
                seq_name = line[0].split(' ')[0]
                newlines.append((seq_name, line[1], line[2], '.', \
                line[4], line[5]))
                ## save as bed6 format later
        if len(newlines) > 0:
            with open(base_dir+ G4_clean_dir + i, 'w') as f0:
                ## substitude GCF with GCA to match GFF files
                f0.write('\n'.join('{}\t{}\t{}\t{}\t{}\t{}'.format(\
                x[0], x[1], x[2], x[3], x[4], x[5]) for x in newlines))
        else:
            continue

def trim_g4_chr_with_seq(base_dir):
    """trim off the extra text to match chromsome name in
    the GFF files to use bedfiles,
    save the new file to the folder all_G4_clean"""
    #base_dir='/Users/Yun/Documents/bacteria_G4/D_thermus/'
    G4_dir = base_dir + "all_G4/"
    if not os.path.isdir(base_dir + 'all_G4_with_seq'):
        os.mkdir(base_dir + 'all_G4_with_seq/')
    for i in os.listdir(G4_dir):
        if i.startswith('.'):
            continue ## ignore the hidden files from apple
        with open(G4_dir+i, 'r') as fp:
            lines = fp.readlines()
            newlines = []
            for line in lines:
                line = line.split('\t')
                seq_name = line[0].split(' ')[0]
                newlines.append((seq_name, line[1], line[2], line[6].split()[0], \
                line[4], line[5]))
                ## save as bed6 format later
        if len(newlines) > 0:
            with open(base_dir+'all_G4_with_seq/' + i, 'w') as f0:
                ## substitude GCF with GCA to match GFF files
                f0.write('\n'.join('{}\t{}\t{}\t{}\t{}\t{}'.format(\
                x[0], x[1], x[2], x[3], x[4], x[5]) for x in newlines))
        else:
            continue

def fasta_parser(fasta_path):
    """parse a fasta file, keep all the seq name as the key
    and the sequence as the value in the seq_dict,
    return the dictionary"""
    g_dict={}
    with open(fasta_path, 'r') as f0:
        lines = f0.readlines()
        seq=''
        for line in lines:
            if line.startswith('>'):
                try: # assign the previous sequence to the previous key
                    g_dict[key]=seq
                except: # handle the first key, not assigned error
                    pass
                key=line.strip()[1:]
                seq=''
            else:
                seq+=line.strip()
    g_dict[key]=seq ## assign the last sequence to the key
    return g_dict
