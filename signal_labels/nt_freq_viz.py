#!/usr/bin/env python 
"""
Map the nucleotide frequency at individual position of a label sequence.

Usage: python nt_freq_viz.py label.fa 
"""

import sys 
from Bio import SeqIO 
from collections import defaultdict


def __main__():
    """
    """

    try:
        fas_file = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1) 

    label_seq_length = 2400  # please change according to the label length 
    nuc_freq_table = defaultdict() 

    for xq in range(label_seq_length):
        nuc_freq_table[xq] = dict(A=0, T=0, C=0, G=0) 
       
    for rec in SeqIO.parse(fas_file, "fasta"):
        for xp, ntd in enumerate(rec.seq):
            nuc_freq_table[xp][ntd] += 1
        
    pos_pref = defaultdict() # max frequency nucleotide in the label sequence positions 
    for xp, freq in nuc_freq_table.items():
        freq_nt = max(freq, key=freq.get) 
        
        pos_pref[xp] = {freq_nt : freq[freq_nt]}
        break  

    # TODO plot function 



if __name__ == "__main__":
    __main__() 
