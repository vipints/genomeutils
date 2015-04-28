#!/usr/bin/env python 
"""
From two different data set generate a set of training and testing labels without the overlap.

Usage: python uniq_train_test_data.py 

Requirement:
    biopython:- http://biopython.org
"""

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

    train_recs = fasta_reader(label_type_1)
    test_recs = fasta_reader(label_type_2)

    shared_keys = set(train_recs.keys()) & set(test_recs.keys())

    print len(shared_keys)


def fasta_reader(fasname):
    """
    read fasta file and return the records locations
    """
    location_marks = dict() 
    for rec in SeqIO.parse(fasname, "fasta"):

        location_marks[rec.id] = 0 

    return location_marks 


if __name__=="__main__":
    __main__()
