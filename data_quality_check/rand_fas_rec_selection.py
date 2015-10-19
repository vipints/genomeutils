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


def fetch_random_example_seq(fasta_in, fa_out, category="+1", total_record_count=46000, sub_sample_records=4000):
    """
    fetch random example signal sequence record from a fasta file using the label type
    
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


def matching_pos_neg_example(in_fas_pos, in_fas_neg):
    """
    fetch random pos examples and corresponding neg example in 1:3 ratio 
    """
    
    from collections import defaultdict 

    total_record_count = 4200 ## based on the pos examples 
    sub_sample_records = 1000  
    neg_sample_ratio = 3 ## sub sample ratio pos:neg  

    try:
        accept_prob = (1.0*sub_sample_records)/total_record_count
    except:
        accept_prob = 1

    print accept_prob
    pos_fas_out = open("tss_sig_pos_example.fa", "w") 
    
    fasta_rec = 0 
    sample_cnt = 1
    pos_fas_rec = defaultdict(list) 
    for rec in SeqIO.parse(in_fas_pos, 'fasta'):
        fasta_rec += 1 
        desc = rec.description.split(" ")
        rnb = random.random()

        if rnb <= accept_prob:
            pos_fas_out.write(rec.format("fasta"))
            pos_fas_rec[desc[-1]].append(rec.description)
            if sample_cnt == sub_sample_records:
                break 
            sample_cnt += 1 

    pos_fas_out.close()

    sys.stdout.write('%d number of records scanned\n' % fasta_rec)
    sys.stdout.write('%d number of examples dumped\n' % sample_cnt)

    pos_neg_ratio = defaultdict(int) 
    neg_fas_out = open("tss_sig_neg_example.fa", "w") 
    for rec_neg in SeqIO.parse(in_fas_neg, 'fasta'):
        desc_neg = rec_neg.description.split(" ")

        if desc_neg[-1] in pos_fas_rec:
            pos_neg_ratio[desc_neg[-1]] += 1
             
            if pos_neg_ratio[desc_neg[-1]] > neg_sample_ratio:
                continue 

            neg_fas_out.write(rec_neg.format("fasta"))

        #import ipdb 
        #ipdb.set_trace()

    neg_fas_out.close() 


if __name__=="__main__": 

    """
    try:
        inh = open(sys.argv[1], "rU")
        outh = open(sys.argv[2], 'w')
    except:
        print __doc__
        sys.exit(-1)

    fetch_random_example_seq(inh, outh)
    """
    
    try:
        in_fasta_pos = sys.argv[1]
        in_fasta_neg = sys.argv[2]
    except:
        exit(__doc__)

    matching_pos_neg_example(in_fasta_pos, in_fasta_neg)

