#!/Users/Yun/anaconda2/bin/python
from collections import Counter
import sys
import string
import argparse ## use for reading file
import os

parser = argparse.ArgumentParser(description = """to filter GC centent of cruciform sequences""", \
                          formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--myfile', '-f',
                   type= str,
                   help='''Input file in fasta format containing one or more
sequences. Use '-' to get the name of the file from stdin

                   ''',
                   required= True)
parser.add_argument('--basedir', '-bd',
                   type= str,
                   help='''The directory to store the output file
                   ''',
                   required= True)

args = parser.parse_args()
## function
def gc(mySeq):
        """calcuate GC percentage of a sequence
        return a float number"""
        mySeq = mySeq.upper()
        return float((mySeq.count("G")+mySeq.count("C"))/float(len(mySeq)))

## now open the file with cruciform DNa and separate the sequence by the GC content
# :
dnaGRich = []
dnaGPoor = []

if args.myfile == '-':
    fh= sys.stdin
else:
    fh= open(args.myfile)
gcLine = 0.41 ## set the cut off GC fraction
lines = fh.readlines()
for line in lines:
    line = line.split()
    seq= line[5]
    gc_percent = gc(seq)
    line.append(gc_percent)
    if gc_percent > gcLine:
        dnaGRich.append(line)
    else:
        dnaGPoor.append(line)
fh.close()
#print(dnaGPoor)
# write into files
#with open( os.path.join(args.basedir, 'GCrich.txt'), 'w') as f0:
    #f0.write('\n'.join('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(\
                    #x[0], x[1], x[2], x[3], x[4], x[5], x[6]) for x in dnaGRich))
#with open( os.path.join(args.basedir, 'GCpoor.txt'), 'w') as f0:
    #f0.write('\n'.join('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(\
                    #x[0], x[1], x[2], x[3], x[4], x[5], x[6]) for x in dnaGPoor))




