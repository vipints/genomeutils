#!/usr/bin/env python
"""
General information about contigs lengths in a FASTA file.

USAGE: python fasta_scan.py in.fasta
"""

import re, sys 
from Bio import SeqIO
from operator import itemgetter
from gfftools import helper 

def __main__():

    try:
        fa_name = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    seq_info = dict()
    fah = helper._open_file(fa_name)

    for rec in SeqIO.parse(fah, "fasta"):
        seq_info[rec.id] = len(rec.seq)
        print rec.id, len(rec.seq)
    fah.close()
    
    print 
    print 'Number of FASTA entries: ', len(seq_info)
    for long_one in sorted(seq_info.items(), key=itemgetter(1), reverse=True):
        print 'Long contig length (bp): ', long_one[0], long_one[1]
        break
    for short_one in sorted(seq_info.items(), key=itemgetter(1)):
        print 'Short contig length (bp): ', short_one[0], short_one[1]
        break
    flength = 0 
    for ele in sorted(seq_info.items(), key=itemgetter(1)):
        flength += ele[1]
    print 'Average length of FASTA contig (bp): ', (flength/len(seq_info))
    print 

if __name__=="__main__":__main__()
