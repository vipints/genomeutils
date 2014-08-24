#!/usr/bin/env python 
"""
Convert GIO format to FASTA format.

Usage: gio_to_fasta.py gio_file_path fasta_file

Requirement:
    biopython : http://biopython.org
"""

import sys, re, os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def __main__():

    try:
        contig_files = os.listdir(sys.argv[1] + '/genome/')
        fasta_file = open(sys.argv[2], "w")
    except:
        print __doc__
        sys.exit(-1)

    if len(contig_files) == 0:
        sys.stderr.write("Empty GIO file found\n")
        sys.exit(-1)

    records = []
    contig = {}
    for chr in contig_files:
        num = re.search(r'Contig(.+)\.flat', chr).group(1)
        contig[int(num)] = chr
    contig_files = []  
    for file in contig:
        ch = open(sys.argv[1] + '/genome/' + contig[file], "rU")
        seq = ''
        for l in ch:
            l = l.strip()
            seq += l
        ch.close()
        chrid = re.search(r'(.+)\.flat', contig[file]).group(1)
        rec = SeqRecord(Seq(seq), id=chrid, description="Created by GIO-to-FASTA tool in mGene.web module")
        records.append(rec)
    SeqIO.write(records, fasta_file, "fasta")
    fasta_file.close()
    records = []

if __name__=="__main__": 
    __main__()
