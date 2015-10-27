import random
import shutil
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict


def split_data(signal, in_fas_pos, sub_sample_records, pos_fas_rec, accept_prob):
    """
    """

    out_fas_pos = "%s_sig_pos_example.bkp" % signal
    pos_fas_out = open(out_fas_pos, "w") 

    sample_cnt = 1
    for rec in SeqIO.parse(in_fas_pos, 'fasta'):
        if rec.id in pos_fas_rec:
            continue

        rnb = random.random()
        if rnb <= accept_prob:
            pos_fas_out.write(rec.format("fasta"))
            pos_fas_rec[rec.id] = 0
            if sample_cnt == sub_sample_records:
                break 
            sample_cnt += 1 

    pos_fas_out.close()

    return sample_cnt, pos_fas_rec


def split_data_rest(signal, in_fas_pos, sub_sample_records, pos_fas_rec):
    """
    """
    out_fas_pos = "%s_sig_pos_example.bkp" % signal
    pos_fas_out = open(out_fas_pos, "w") 

    sample_cnt = 1
    for rec in SeqIO.parse(in_fas_pos, 'fasta'):
        if rec.id in pos_fas_rec:
            continue

        pos_fas_out.write(rec.format("fasta"))
        pos_fas_rec[rec.id] = 0
        if sample_cnt == sub_sample_records:
            break 
        sample_cnt += 1 

    pos_fas_out.close()


def split_data_random_non_overlap(signal="tss", total_record_count=5000):
    """
    """

    in_fas_pos = "%s_sig_pos_example.fa" % signal
    #FIXME 
    sub_sample_records = 1000

    try:
        accept_prob = (1.0*sub_sample_records)/total_record_count
    except:
        accept_prob = 1
    
    while True:
        pos_fas_rec = defaultdict()
        counter, pos_fas_rec = split_data(signal, in_fas_pos, sub_sample_records, pos_fas_rec, accept_prob)

        if counter == sub_sample_records:
            break 
        
        print "set_1 recounting" 
    shutil.move('%s_sig_pos_example.bkp' % signal, '%s_sig_pos_example_1.fa' % signal)

    total_record_count -= sub_sample_records
    try:
        accept_prob = (1.0*sub_sample_records)/total_record_count
    except:
        accept_prob = 1

    accept_prob += 0.1

    while True:
        counter, pos_fas_rec = split_data(signal, in_fas_pos, sub_sample_records, pos_fas_rec, accept_prob)

        if counter == sub_sample_records:
            break 

        print "set_2 recounting" 
    shutil.move('%s_sig_pos_example.bkp' % signal, '%s_sig_pos_example_2.fa' % signal)

    #FIXME 
    sub_sample_records = 3000
    split_data_rest(signal, in_fas_pos, sub_sample_records, pos_fas_rec)

    shutil.move('%s_sig_pos_example.bkp' % signal, '%s_sig_pos_example_3.fa' % signal)

    #FIXME
    # need a master module to control the cleaning of pos and neg example cleanup 


def get_matching_neg_example(in_fas_pos, in_fas_neg, out_fas_neg):
    """
    """

    #in_fas_pos = "%s_sig_pos_example.fa" % signal
    #in_fas_neg = "%s_sig_neg_example.fa" % signal
    
    pos_fas_rec = defaultdict(list) 
    for rec in SeqIO.parse(in_fas_pos, 'fasta'):
        desc = rec.description.split(" ")

        pos_fas_out.write(rec.format("fasta"))
        pos_fas_rec[desc[-1]].append(rec.description) 
        
    print "%d number of positive records" % len(pos_fas_rec) 

    neg_sample_ratio = 3 ## sub sample ratio pos:neg  
    neg_rec_cnt = 0 
    pos_neg_ratio = defaultdict(int) 

    neg_fas_out = open(out_fas_neg, "w") 
    for rec_neg in SeqIO.parse(in_fas_neg, 'fasta'):
        desc_neg = rec_neg.description.split(" ")

        if desc_neg[-1] in pos_fas_rec:
            pos_neg_ratio[desc_neg[-1]] += 1
             
            if pos_neg_ratio[desc_neg[-1]] > neg_sample_ratio:
                continue 

            neg_fas_out.write(rec_neg.format("fasta"))
            neg_rec_cnt += 1 

    print "%d number of negative records" % neg_rec_cnt 
    neg_fas_out.close()

