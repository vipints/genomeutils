#!/usr/bin/env python
"""
Program to create a stable release of genome sequence in FASTA format. 
It filters the partially assembled contigs according to the users choice.

Usage:
cat contigs.txt | python make_stable_genome.py in.fasta out.fasta 

contigs.txt file looks like:
Chr1
Chr2
Chr3
Chr4

Please adjust your valid contigs file to create a subset valid file. 

Requirement:
    BioPython:- http://biopython.org 
    gfftools:- https://github.com/vipints/genomeutils/tree/master/gfftools/helper.py 
"""
import sys, os
from Bio import SeqIO 
from gfftools import helper

def stop_err(msg):
    """
    stop the execution and print out the captured error message. 
    """
    sys.stderr.write('%s\n' % msg)
    sys.exit(-1)

def __main__():
    
    try:
        fname = sys.argv[1]
        fa_out = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)

    # get the valid chromosome identifier from user as STDIN
    chrom = dict()
    for chr in sys.stdin:
        chr = chr.strip()
        chrom[chr] = 1
    
    # get the filehandler from input file
    try:
        fh = helper._open_file(fname)
    except Exception, errmsg:
        stop_err('error in reading file '+ errmsg) 

    # check the out filehandler
    try:
        outfh = open(fa_out, "w")
    except Exception, errmsg:
        stop_err('error in writing file '+ errmsg) 

    # writing stable contig genome sequence in FASTA format 
    for rec in SeqIO.parse(fh, "fasta"):
        if rec.id in chrom:
            outfh.write(rec.format("fasta"))
    fh.close()
    outfh.close()

if __name__=="__main__":
    __main__()
