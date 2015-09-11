#!/usr/bin/env python
"""
select training/testing dataset labels according to positive or negative class 
from ARTS gold standard dataset for transcript start sites.

Usage:
    rand_fas_rec_selection.py in_fasta.fa out_fasta.fa

Requirement:
    biopython: http://biopython.org 
"""

import sys
import random
from Bio import SeqIO


def fetch_random_label_seq(fasta_in, fa_out, category="+1", total_record_count=46000, sub_sample_records=4000):
    """
    fetch random label sequence record from a fasta file using the label type
    
    @args fasta_in: fastafile with label records 
    @type fasta_in: str or a filehandler 
    @args fa_out: outfile handler 
    @type fa_out: str 
    @args category: label class 
    @type category: str 
    """

    label_type = ["-1", "+1"]
    if not category in label_type:
        sys.exit("error: sequence record label type +1/-1 not supported to %s\n" % category)

    try:
        accept_prob = (1.0*sub_sample_records)/total_record_count
    except:
        accept_prob = 1
    
    accept_prob = 0.939 
    cnt = 1
    fas_rec = 0 

    for rec in SeqIO.parse(fasta_in, 'fasta'):
        desc = rec.description.split(" ")

        if desc[1]==category:
            rnb = random.random()

            if rnb <= accept_prob:
                fa_out.write(rec.format("fasta"))

                if cnt == sub_sample_records:
                    break 

                cnt += 1 
        fas_rec += 1 

    fasta_in.close()
    fa_out.close()

    sys.stdout.write('\n%d number of records scanned\n' % fas_rec)
    sys.stdout.write('%d number of labels dumped\n\n' % cnt)


if __name__=="__main__": 

    try:
        inh = open(sys.argv[1], "rU")
        outh = open(sys.argv[2], 'w')
    except:
        print __doc__
        sys.exit(-1)
    
    fetch_random_label_seq(inh, outh)
