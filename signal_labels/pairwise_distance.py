#!/usr/bin/env python 
"""
Computing the distance between pairs of sequences which are derived from a multiple sequence alignment. 

    TODO 
Requirement:
    biopython :- http://biopython.org/
    scikit-bio :- https://github.com/biocore/scikit-bio
"""


from Bio import AlignIO, SeqIO
from skbio import RNA, Alignment 


def compute_distance_matrix(msa_file):
    """
    load up some aligned sequences, and compute a distance matrix 
    compute distances between the sequences using the hamming function

    see also: 
    scipy.spatial.distance.hamming

    @args msa_file: multiple sequence alignment in fasta format 
    @type msa_file: str 

    """
    
    records = [] 
    for rec in SeqIO.parse(msa_file, "fasta"):
        records.append(RNA(rec.seq, rec.id))

    aln = Alignment(records)
    master_dm = aln.distances() 

    ## writing the result to a csv file 
    csvfile = "distance_mat.csv"

    ## result as a list of list 
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerows(master_dm)

