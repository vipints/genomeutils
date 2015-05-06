#!/usr/bin/env python 
"""
Computing the distance between pairs of sequences which are derived from a multiple sequence alignment. 

    
Requirement:
    biopython :- http://biopython.org/
    scikit-bio :- https://github.com/biocore/scikit-bio
    NINJA :- http://nimbletwist.com/software/ninja 
"""

import os 
import csv 
import subprocess
from Bio import AlignIO, SeqIO
from skbio import RNA, Alignment 


def compute_distance_matrix(msa_file, csvfile="distance_mat.csv"):
    """
    load up some aligned sequences, and compute a distance matrix 
    compute distances between the sequences using the hamming function

    see also: 
    scipy.spatial.distance.hamming

    @args msa_file: multiple sequence alignment in fasta format 
    @type msa_file: str 
    @args csvfile: output distance matrix file in csv format 
    @type csvfile: str 
    """
    
    records = [] 
    for rec in SeqIO.parse(msa_file, "fasta"):
        records.append(RNA(rec.seq, rec.id))

    aln = Alignment(records)
    master_dm = aln.distances() 

    ## writing the result to a csv file 
    csv_header_row = [header for header in master_dm.ids] 

    ## result as a list of list 
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator="\n")

        writer.writerows([csv_header_row])
        writer.writerows(master_dm)
    
    output.close() 


def run_ninja(msa_file, distance_mat="distance_matrix.csv"):
    """
    NINJA is software for inferring large-scale neighbor-joining phylogenies.

    @args msa_file: multiple sequence alignment in fasta format 
    @type msa_file: str 
    @args distance_mat: output distance matrix file in csv format 
    @type distance_mat: str 

    NOTE: This module expects ninja is available under PATH variable or add the path below line.
    """

    #os.environ['PATH'] += os.pathsep + '/home/share/software/ninja/'

    cli = 'ninja --alph_type d --out_type d --corr_type n %s > %s' % (msa_file, distance_mat) 
    process = subprocess.Popen(cli, shell=True) 
    process.wait()
