#!/Users/Yun/anaconda3/bin/python
import re
import sys
import string
import argparse ## adding parameters
import operator ## use for the sorting tabl
import multiprocessing as mp ## use for parallel programming

## take two arguments, the fasta file and the number of processors to use
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
    """take a sequence and return the reverse complementary strand of this seq
    input: a DNA sequence
    output: the complementary of this seq"""
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'O'}
    ##get the complementary if key is in the dict, otherwise return the base, set O base pair with N to remove NNNN in the genome
    return(''.join(bases_dict.get(base, base) for base in seq[::-1]))
print([x for x in range(12,5,-1)])


def get_cruci_list(seq_name, seq):
    """take a sequence and return the location, loop size, stem size of the cruciform DNA
    todo: vectorize with numpy"""
    cruci = []
    pos = len(seq)-(5+12+1) ## the last position to check
    max_loop =8
    max_stem = 11
    t = max_stem + max_loop ## starting point, t sit in the middle of the cruciform structure
    while t < pos:
        jump=False
        for i in range(3,8,1): ## the even number loop
            l=int(i/2)
            for j in range(11, 5, -1): ## for stem 18 to 6
                if i%2==0:
                    if seq[(t-j-l):(t-l)].upper() == comp(seq[(t+l):(t+l+j)].upper()):
                        cruci.append([seq_name, t-j-l,t+l+j, i, j, seq[t-j-l:t+j+l]])
                    #cruci_set.update([t+j, i]) ## add the start and the end position to the set, make sure it is unique
                        t = t + 2*j + i ## t=t+j+4+6 assume the small length?
                        jump=True
                        break
                elif i%2!=0: #odd number loop
                    if seq[(t-j-l):(t-l)].upper() == comp(seq[(t+l+1):(t+l+1+j)].upper()):
                        cruci.append([seq_name, t-j-l,t+l+j+1, i, j, seq[t-j-l:t+j+l+1]]) ## add one extra base for the l=int(i/2) step
                    #cruci_set.update([t+j, i]) ## add the start and the end position to the set, make sure it is unique
                        t = t + 2*j + i ## t=t+j+4+6 assume the small length?
                        jump=True
                        break
            if jump: ## find a match, stop looking for the smaller one
                    break # break from the loop scan

            elif i==4: ## the last loop length, add 1 to cursor
                t+=1 # move the cursor 1 bp if didn't find any match

    return(cruci)

def fasta_parser(f0):
    """parse a fasta file, keep all the seq name as the key
    and the sequence as the value in the seq_dict,
    return the dictionary
    input: f0 is the file handle
    output: g_dict is the sequence dictionary"""
    g_dict={}
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

##-------------------------------------read and process files line by line-----------------
if args.fasta == '-':
    ref_seq_fh= sys.stdin
else:
    ref_seq_fh= open(args.fasta)
cruci_all = []
ref_seq=fasta_parser(ref_seq_fh)
for key, seq in ref_seq.items():
    cruci_all.extend(get_cruci_list(key,seq)) ## extend the list of cruci to all cruci

ref_seq_fh.close() ## close the handle when done

for line in cruci_all:
    line= '\t'.join([str(x) for x in line])
    print(line)

sys.exit()