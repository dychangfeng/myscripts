#!/Users/Yun/anaconda2/bin/python
import re
import sys
import string
import argparse ## adding parameters
import operator ## use for the sorting tabl
from datetime import datetime

start_at=datetime.now()
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

ref_seq_fh.close() ## close the handle when done

#cruci_final= sorted(cruci_all, key=operator.itemgetter(0,1,2)) # sort by name, start and end position

for line in cruci_all:
    line= '\t'.join([str(x) for x in line])
    print(line)

total_time=datetime.now()-start_at
print(total_time)

sys.exit()