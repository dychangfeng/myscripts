#!/Users/Yun/anaconda2/bin/python
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
parser.add_argument('--processors', '-p',
                   type= int,
                   help='''the number of processors to used for this calculation

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


def get_cruci_list(seq_name, seq, i):
    """take a sequence and return the location, loop size, stem size of the cruciform DNA
    i is the loop length that examed
    todo: vectorize with numpy"""
    cruci = []
    pos = len(seq)-(5+12+1) ## the last position to check
    t = 0
    while t < pos:
        for j in range(12, 5, -1): ## for stem 15 to 6
            if seq[t:(t+j)].upper() == comp(seq[(t+j+i):(t+2*j+i)].upper()):
                cruci.append([seq_name, t,t+2*j+i, i, j, seq[t:t+2*j+i]])
                #cruci_set.update([t+j, i]) ## add the start and the end position to the set, make sure it is unique
                t = t + 2*j + i

            elif j<6: ## until the last loop, no match
                t+=1 # move the cursor 1 bp
    return(cruci)

def multi_loops(n_processes, seq_name, seq, loop_lens):
    """n_processes is the number of precessors to use
    seq_name, seq is the name and sequence of the sequence the be checked
    loop_lens is a list of loop lengths to be examed"""
    pool = mp.Pool(processes=n_processes)
    results = [pool.apply_async(get_cruci_list, args=(seq_name, seq, i)) for i in loop_lens]
    final = []
    for p in results:
        final.extend(p.get()) ## merge to one list
    final_sorted= sorted(final, key=operator.itemgetter(0,1,2)) # sort by name, start and end position
    return final_sorted
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
loops = [6,7] ## define the loops of the cruciform DNA
while True:
    while line.startswith('>') is False:
        ## get the list of ref_seq for this chr
        ref_seq.append(line)
        line= (ref_seq_fh.readline()).strip()
        if line == '':
            break
    ref_seq= ''.join(ref_seq) # join the sequence of the ref seq
    cruci_all.extend(multi_loops(args.processors, chr, ref_seq,loops)) ## extend the list of cruci to all cruci
    chr= re.sub('^>', '', line)
    ref_seq= [] ## clear the ref_seq
    line= (ref_seq_fh.readline()).strip() ## read a new line that it does not start with a ">"
    if line == '':
        break

ref_seq_fh.close() ## close the handle when done

for line in cruci_all:
    line= '\t'.join([str(x) for x in line])
    print(line)

sys.exit()