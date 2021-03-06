#!/usr/bin/env python
"""
Visual check to the label sequences used for training and testing SVMs.

Usage: python seq_consensus_viz.py in.fasta fig1.pdf 

Requirements:
    BioPython:- http://biopython.org 
    pylab:- http://wiki.scipy.org/PyLab 
    numpy:- http://numpy.org 
"""

import sys
from Bio import SeqIO


def nuc_freq_viz(fasta_file, label_seq_length=2401):
    """
    nucleotide frequency visualization 

    @args fasta_file: label sequence in fasta format 
    @type fasta_file: str 
    @args label_seq_length: number of nucleotide in a label sequence (default 2401)
    @type label_seq_length: int
    """

    from collections import defaultdict

    nuc_freq_table = defaultdict() 

    for xq in range(int(label_seq_length)):
        nuc_freq_table[xq] = dict(A=0, T=0, C=0, G=0) 
       
    for rec in SeqIO.parse(fasta_file, "fasta"):
        for xp, ntd in enumerate(rec.seq):
            nuc_freq_table[xp][ntd] += 1

    pos_pref = defaultdict() # max frequency nucleotide in the label sequence positions 

    for xp, freq in nuc_freq_table.items():

        freq_nt = max(freq, key=freq.get) 
        pos_pref[xp] = {freq_nt : freq[freq_nt]}
        break  
        #TODO 

    # TODO plot function 


def seq_viz_stack(fasta_file, pdf_file_name):
    """
    color coded sequence visualization 

    @args fasta_file: label sequences in fasta format
    @type fasta_file: str 
    @args pdf_file_name: plotting result file 
    @type pdf_file_name: str 
    """

    seq_records = [] 
    for rec in SeqIO.parse(fasta_file, "fasta"):
        if rec.seq:
            seq_records.append(str(rec.seq))

    import numpy 
    import pylab as pl 

    pl.figure(figsize=(30,len(seq_records)/100))

    num_seqs = len(seq_records)
    len_seqs = len(seq_records[0])

    color_matrix = numpy.zeros((num_seqs, len_seqs))

    base_map = dict(A = 0, T = 1, G = 2, C = 3, N = 4) 

    for ix, sequence in enumerate(seq_records):
        for ij, nucleotide in enumerate(sequence):
            color_matrix[ix,ij] = base_map[nucleotide]

    pl.imshow(color_matrix)

    #TODO add the legend to the plot

    pl.title("sequence distribution")
    pl.xlabel("sequence position")
    pl.ylabel("sequence")

    pl.savefig(pdf_file_name) 


if __name__ == "__main__":

    try:
        fas_name = sys.argv[1]
        res_name = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)

    seq_viz_stack(fas_name, res_name) 
