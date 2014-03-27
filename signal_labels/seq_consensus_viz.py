#!/usr/bin/env python
"""
Visual check to the label sequences.

Usage: python seq_consensus_viz.py in.fasta

Requirements:
    BioPython:- http://biopython.org 
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def seq_viz(records):
    """
    color coded sequence visualization 
    """

    import numpy 
    import pylab as pl 

    pl.figure()

    num_seqs = len(records)
    len_seqs = len(records[0])

    color_matrix = numpy.zeros((num_seqs, len_seqs))

    base_map = dict(A = 0, T = 1, G = 2, C = 3) 

    for ix, sequence in enumerate(records):
        for ij, nucleotide in enumerate(sequence):
            color_matrix[ix,ij] = base_map[nucleotide]

    pl.imshow(color_matrix)

    pl.title("sequence distribution")
    pl.xlabel("sequence position")
    pl.ylabel("sequence")

    pl.show() 
    

def __main__():

    try:
        fas_name = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    
    seq_records = [] 
    for rec in SeqIO.parse(fas_name, "fasta"):
        if rec.seq:
            seq_records.append(str(rec.seq))
    
    seq_viz(seq_records) 


if __name__ == "__main__":
    __main__()
