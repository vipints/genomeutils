#!/usr/bin/env python 
"""Find the sequence similarity ratio of two sequences. 
"""
import re, sys 
from Bio import SeqIO 
import difflib

def read_fasta(fname):
    frec=dict()
    fah=open(fname, 'rU')
    for rec in SeqIO.parse(fah, 'fasta'):
        frec[rec.id]=str(rec.seq) 
    fah.close()
    return frec

try:
    fas_1=sys.argv[1]
    fas_2=sys.argv[2]
except:
    print __doc__
    sys.exit(-1)
    
ele_1=read_fasta(fas_1)
ele_2=read_fasta(fas_2)

for fid, fseq in sorted(ele_1.items()):
    for fsd, eachseq in sorted(ele_2.items()):
        eachseq=re.sub(r'U', r'T', eachseq)
        if difflib.SequenceMatcher(None, fseq, eachseq).ratio()==1.0:
            print fid, fsd
        
