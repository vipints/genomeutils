#!/usr/bin/env python 
"""
From two different data set generate a set of training and testing labels without the overlap.

Usage: python uniq_train_test_data.py 

Requirement:

"""

import re 
import sys
from Bio import SeqIO 

def __main__():
    """
    """

    try:
        label_type_1 = sys.argv[1] # train data 
        label_type_2 = sys.argv[2] # test data 
    except: 
        print __doc__
        sys.exit(-1) 


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
