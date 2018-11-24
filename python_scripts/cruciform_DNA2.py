#!/Users/Yun/anaconda2/bin/python
import re
import sys
import string
import argparse ## use for reading file
import operator ## use for the sorting tabl
import multiprocessing as mp ## use for parallel programming

parser = argparse.ArgumentParser(description = """to find cruciform DNA""",formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-f',
                   type= str,
                   help='''Input file in fasta format containing one or more
sequences. Use '-' to get the name of the file from stdin

                   ''',
                   required= True)


args = parser.parse_args()
##-------------------------functions-------------------------
## inverted repeats are reverse complementary of each other
def comp(seq):
    """take a sequence and return the complementary strand of this seq
    input: a DNA sequence
    output: the complementary of this seq"""
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'O'}
    ##get the complementary if key is in the dict, otherwise return the base, set O base pair with N to remove NNNN in the genome
    return(''.join(bases_dict.get(base, base) for base in seq[::-1]))



def get_cruci_list(seq_name, seq):
    """take a sequence and return the location, loop size, stem size of the cruciform DNA
    todo: vectorize with numpy"""
    cruci = []
    cruci_set = set() ## store the unique value of the loop and start position to avoid duplicates
    pos = len(seq)-(5+12+1) ## the last position to check
    for i in range(4,6): # loop size 4, 5
        print('sceening loop size {} for sequence {}'.format(i, seq_name))
        t = 0
        while t < pos:
            for j in range(12, 5, -1): ## for stem 12 to 8
                if seq[t:(t+j)].upper() == comp(seq[(t+j+i):(t+2*j+i)].upper()):
                    cruci.append([seq_name, t,t+2*j+i, i, j, seq[t:t+2*j+i]])
                        #cruci_set.update([t+j, i]) ## add the start and the end position to the set, make sure it is unique
                    t = t + 2*j + i # move cursor to the end of the hairpin， scan window
                else: ## until the last loop, no match
                    t+=1 # move the cursor 1 bp
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