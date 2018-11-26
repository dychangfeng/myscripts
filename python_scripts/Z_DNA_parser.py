#!/usr/bin/python

import re
import sys
import string
import argparse
import operator

parser = argparse.ArgumentParser(description="""

DESCRIPTION

    Search for matches to a regex in a fasta file and return a bed file with the
    coordinates of the match and the matched sequence itself.

    With defaults, quadparser.py searches for putative quadruplexes on forward
    and reverse strand using the quadruplex rule described at
    http://en.wikipedia.org/wiki/G-quadruplex#Quadruplex_prediction_techniques.

    The defualt regex is '([Cc][gG]){5,}' and its complement
    produce the same output as in http://www.quadruplex.org/?view=quadbaseDownload

    Output bed file has columns:
    1. Name of fasta sequence (e.g. chromosome)
    2. Start of the match
    3. End of the match
    4. ID of the match
    5. Length of the match
    6. Strand
    7. Matched sequence

    Note: Fasta sequences are read in memory one at a time. Also the bed file
    of of the matches are kept in memeory.

    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--regex', '-r',
                   type= str,
                   help='''Regex to be searched in the fasta input.
Matches to this regex will have + strand. This string passed to python
re.compile(). The default regex is '([Cc][gG]){5,}' which searches
for Z-DNA.

                   ''',
                   default= '([Cc][gG]){5,}')

parser.add_argument('--regexrev', '-R',
                   type= str,
                   help='''The second regex to be searched in fasta input.
Matches to this regex will have - strand.
By default (None), --regexrev will be --regex complemented by replacing
'actguACTGU' with 'tgacaTGACA'. NB: This means that [a-zA-Z] will be translated
to [t-zT-Z] and proteins are not correctly translated.

                   ''',
                   default= None)


parser.add_argument('--fasta', '-f',
                   type= str,
                   help='''Input file in fasta format containing one or more
sequences. Use '-' to get the name of the file from stdin

                   ''',
                   required= True)

parser.add_argument('--noreverse',
                   action= 'store_true',
                   help='''Do not search the reverse (-) strand. I.e. do not use
the complemented regex (or --regexrev/-R). Use this flag to search protein
sequences.

                   ''')

args = parser.parse_args()

" --------------------------[ Check and pare arguments ]---------------------- "

""" Reverse forward match """
intab=  'actguACTGU'
outtab= 'tgacaTGACA'
if args.regexrev is None:
    transtab = str.maketrans(intab, outtab)
    regexrev= args.regex.translate(transtab)
else:
    regexrev= args.regex
# handle stdin
if args.fasta == '-':
    args.fasta= sys.stdin.readlines()
    if len(args.fasta) > 1:
        sys.exit('\n Only one input file at a time can be processed:\n--fasta/-f: %s\n' %(args.fasta))
    args.fasta= args.fasta[0].strip()

" ------------------------------[  Functions ]--------------------------------- "

def reverse_comp(seq):
    """take a sequence and return the reverse complementary strand of this seq
    input: a DNA sequence
    output: the reverse complementary of this seq"""
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    ##get the complementary if key is in the dict, otherwise return the base
    return ''.join(bases_dict.get(base, base) for base in seq[::-1])

psq_re_f= re.compile(args.regex)
psq_re_r= re.compile(regexrev)

ref_seq_fh= open(args.fasta)

ref_seq=[]
line= (ref_seq_fh.readline()).strip()
chr= re.sub('^>', '', line)
line= (ref_seq_fh.readline()).strip()
gquad_list= []
while True:
    while line.startswith('>') is False:
        ref_seq.append(line)
        line= (ref_seq_fh.readline()).strip()
        if line == '':
            break
    ref_seq= ''.join(ref_seq)
    for m in re.finditer(psq_re_f, ref_seq):
        quad_id= str(chr) + '_' + str(m.start()) + '_' + str(m.end()) + '_for'
        gquad_list.append([chr, m.start(), m.end(), quad_id, len(m.group(0)), '+', reverse_comp(m.group(0).upper())])
    if args.noreverse is False:
        for m in re.finditer(psq_re_r, ref_seq):
            quad_id= str(chr) + '_' + str(m.start()) + '_' + str(m.end()) + '_rev'
            gquad_list.append([chr, m.start(), m.end(), quad_id, len(m.group(0)), '-', reverse_comp(m.group(0).upper())])
    chr= re.sub('^>', '', line)
    ref_seq= []
    line= (ref_seq_fh.readline()).strip()
    if line == '':
        break

gquad_sorted= sorted(gquad_list, key=operator.itemgetter(0,1,2,3))

for line in gquad_sorted:
    line= '\t'.join([str(x) for x in line])
    print(line)
sys.exit()
