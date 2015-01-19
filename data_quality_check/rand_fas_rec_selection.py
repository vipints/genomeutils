#!/usr/bin/env python
"""
select random fasta records from a fasta file. 

"""

import sys
import random
from Bio import SeqIO


def __main__():

    try:
        inh = open(sys.argv[1], "rU")
        fasta_out = open(sys.argv[2], 'w')
    except:
        print __doc__
        sys.exit(-1)

    total_record_count = 46000
    sub_sample_records = 12000 # 12000

    label_type = "-1" # -1 

    try:
        accept_prob = (1.0*sub_sample_records)/total_record_count
    except:
        accept_prob = 1
    
    #accept_prob = 0.285
    print accept_prob

    cnt = 1
    fas_rec = 0 

    for rec in SeqIO.parse(inh, 'fasta'):
        desc = rec.description.split(" ")

        if desc[1]==label_type:
            rnb = random.random()

            if rnb <= accept_prob:
                fasta_out.write(rec.format("fasta"))

                if cnt == sub_sample_records:
                    break 

                cnt += 1 
        fas_rec += 1 

    inh.close()
    fasta_out.close()

    print
    print '%d number of records scanned' % fas_rec 
    print '%d number of labels dumped' % cnt 
    print 

if __name__=="__main__": 
    __main__()
