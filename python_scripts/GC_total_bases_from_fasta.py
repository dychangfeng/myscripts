#!/usr/bin/python

import sys
import string
import argparse

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    calculate GC content of every fasta files and the length of each genome
    """, formatter_class= argparse.RawTextHelpFormatter)
    
parser.add_argument('--fasta', '-f',
                   type= str,
                   help='''Input file in fasta format containing one or more
sequences. Use '-' to get the name of the file from stdin                   
                   ''',
                   required= True)
args = parser.parse_args()

if args.fasta == '-':
    args.fasta= sys.stdin.readlines()
    if len(args.fasta) > 1:
        sys.exit('\nquadpareser.py: Only one input file at a time can be processed:\n--fasta/-f: %s\n' %(args.fasta))
    args.fasta= args.fasta[0].strip()
ref_seq_fh= open(args.fasta)

seq = ''
lst=list()

for line in ref_seq_fh:
# let's discard the newline at the end (if any)
	line=line.rstrip()
# distinguish header from sequence
	if line[0]!='>': # or line.startswith('>')
		line=line.replace(" ", "")
		line=line.replace("\n", "")
		seq = seq + line
		#print seqs[name]

def GC_percentage(seq):	
	count=0
	total=0
	for base in seq:
		if base == 'C' or base =='G':
			count+=1	
		total+=1
	percentage = count/float(total)
	return percentage, total

print GC_percentage(seq)