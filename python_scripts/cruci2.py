#!/Users/Yun/anaconda2/bin/python
import re
import sys
import string
import argparse ## use for reading file
import operator ## use for the sorting table

parser = argparse.ArgumentParser(description = """to find cruciform DNA""",formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-f',
                   type= str,
                   help='''Input file in fasta format containing one or more
sequences. Use '-' to get the name of the file from stdin

                   ''',
                   required= True)


args = parser.parse_args()
##-------------------------functions-------------------------

def comp(seq):
    """take a sequence and return the complementary strand of this seq
    input: a DNA sequence
    output: the complementary of this seq"""
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    ##get the complementary if key is in the dict, otherwise return the base
    return(''.join(bases_dict.get(base, base) for base in seq[::-1]))


def get_cruci_list(seq_name, seq):
    """take a sequence and return the location, loop size, stem size of the cruciform DNA"""
    cruci = []
    cruci_set = set() ## store the unique value of the loop start position to avoid duplicates
    pos = set() 
    for i in range(3,6): # loop size
        print('screening loop size {}'.format(i))
        for j in range(12,7, -1): # stem size
            print('screening stem size {}'.format(j))
            for t in range(len(seq)-(i+j+1)) : ## loop through the whole sequence
                if seq[t:t+j].upper() == comp(seq[t+j+i:t+2*j+i].upper()):
                    if not set([t, t+j, t+2*j+i]).intersection(cruci_set): 
                    ## check whether the start and the end in our list 
                        cruci.append([seq_name, t,t+2*j+i, i, j, seq[t:t+2*j+i]])
                        cruci_set.update([t, t+j, t+2*j+i]) ## add the start and the end position to the set, make sure it is unique
                    
    return(cruci)


##-------------------------------------read and process files line by line-----------------
if args.fasta == '-':
    ref_seq_fh= sys.stdin
else:
    ref_seq_fh= open(args.fasta)
cruci_all = []
ref_seq=[]
line= (ref_seq_fh.readline()).strip()
chr= re.sub('^>', '', line) ## get the chr name, remove the '>' sign
line= (ref_seq_fh.readline()).strip() # read another line
tt = 0
while True:
    while line.startswith('>') is False:
        ## get the list of ref_seq for this chr
        ref_seq.append(line)
        line= (ref_seq_fh.readline()).strip()
        if line == '':
            break
    ref_seq= ''.join(ref_seq) # join the sequence of the ref seq
    cruci_all.extend(get_cruci_list(chr, ref_seq)) ## extend the list of cruci to all cruci
    chr= re.sub('^>', '', line)
    ref_seq= [] ## clear the ref_seq
    line= (ref_seq_fh.readline()).strip() ## read a new line that it does not start with a ">"
    if line == '':
        break
cruci_sorted= sorted(cruci_all, key=operator.itemgetter(0,1,2)) # sort by name, start and end position
for line in cruci_sorted:
    line= '\t'.join([str(x) for x in line])
    print(line)

sys.exit()