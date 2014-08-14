#!/usr/bin/env python 
"""
From two different data set generate a set of training and testing labels without the overlap.

Usage: python uniq_train_test_data.py 

Requirement:

"""

import re 
import sys
from Bio import SeqIO 
from gfftools import GFFParser

def __main__():
    """
    """

    try:
        label_type_1 = sys.argv[1] # train data 
        label_type_2 = sys.argv[2] # test data 
    except: 
        print __doc__
        sys.exit(-1) 

    #gff_fname = "/cbio/grlab/nobackup/projects/SignalPrediction/SRA-rnaseq/H_sapiens/trans_pred/predgenes.gff3" 
    gff_fname = "/cbio/grlab/home/vipin/tmp/test-data/splice_gene.gff"
    gff_content = GFFParser.Parse(gff_fname)

    for rec in gff_content:
        for xp, sub_rec in enumerate(rec['transcripts']):
            print xp, sub_rec, rec['exons'][xp][0]

    #fasta_reader(label_type_1)


def fasta_reader(fasname):
    """
    read fasta file and return the records locations
    """

    for rec in SeqIO.parse(fasname, "fasta"):
        print rec.id

        break 


if __name__=="__main__":
    __main__()
