#!/usr/bin/env python
"""
Program to find the correct number of valid labels for each type of signals

Usage: delete_label_entity.py H_sapiens/ens_98/ [acc|don|tis|tss] 
"""
import os, re, sys
from Bio import SeqIO

def __main__():
    try:
        base_path = sys.argv[1]
        signal = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1)
   
    out_path = os.path.dirname(base_path) 

    fasta_out_plus = open(out_path + "/" + signal +"_sig_plus_label.bkp", 'w')
    cnt = 0 
    
    #plus_label_cnt = 10000
    plus_label_cnt = 1000
    plus_label_cnt = 1000 if signal in ['acc', 'don'] else plus_label_cnt

    plus_hd = open(out_path + "/"+ signal + "_sig_plus_label.fa", "rU")
    for rec in SeqIO.parse(plus_hd, 'fasta'):
        if not rec.seq:continue
        SeqIO.write([rec], fasta_out_plus, "fasta")
        cnt += 1
        if cnt == plus_label_cnt:
            break
    plus_hd.close()
    fasta_out_plus.close()
    
    os.system('mv ' + out_path + '/' + signal +'_sig_plus_label.bkp '+ out_path + "/"+ signal + "_sig_plus_label.fa")

    fasta_out_minus = open(out_path + "/"+ signal + "_sig_minus_label.bkp", 'w')
    cnt = 0 
    
    #minus_label_cnt = 30000
    minus_label_cnt = 3000
    minus_label_cnt = 3000 if signal in ['acc', 'don'] else minus_label_cnt

    minus_hd = open(out_path + "/" + signal +"_sig_minus_label.fa", "rU")
    for rec in SeqIO.parse(minus_hd, 'fasta'):
        if not rec.seq:continue
        SeqIO.write([rec], fasta_out_minus, "fasta")
        cnt += 1
        if cnt == minus_label_cnt:
            break
    minus_hd.close()
    fasta_out_minus.close()
    
    os.system('mv ' + out_path + '/'+ signal + '_sig_minus_label.bkp '+ out_path + "/" + signal + "_sig_minus_label.fa")
    
    #summary of entries
    cnt = 0 
    signal_seq_len = 401 #default TIS TSS label sequence length 
    signal_seq_len = 200 if signal in ['acc', 'don'] else signal_seq_len

    plus_hd = open(out_path + "/" + signal + "_sig_plus_label.fa", "rU")
    for rec in SeqIO.parse(plus_hd, 'fasta'):
        if not rec.seq:continue
        if len(rec.seq)==signal_seq_len:
            cnt += 1
    plus_hd.close()
    print '+', cnt

    cnt = 0 
    minus_hd = open(out_path + "/" + signal + "_sig_minus_label.fa", "rU")
    for rec in SeqIO.parse(minus_hd, 'fasta'):
        if not rec.seq:continue
        if len(rec.seq)==signal_seq_len:
            cnt += 1
    minus_hd.close()
    print '-', cnt

if __name__=="__main__": __main__()
