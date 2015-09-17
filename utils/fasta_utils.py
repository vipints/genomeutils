#!/usr/bin/env python 

from Bio import SeqIO 
from operator import itemgetter

from gfftools import helper 


def split_fasta_records(fas_fname):
    """
    split the fasta records to multiple fasta file 
    """
    
    cnt = 0 
    records = 100 
    out_file_prefix = "tss_sig_pos_example_"

    for rec in SeqIO.parse(fas_fname, "fasta"):
     
        if cnt % records == 0:
            fout_name = "%s%d.fa" % (out_file_prefix, cnt) 

            try:
                fh.close()
            except:
                pass 
            try: 
                fh = open(fout_name, "w")
            except OSError:
                exit("error: cannot create file %s " % fout_name)

            fh.write(rec.format("fasta"))
        else:
            fh.write(rec.format("fasta"))
                    
        cnt +=1 


def fasta_seq_length(fa_name):
    """
    general information about contigs lengths in a FASTA file
    """

    seq_info = dict()
    fah = helper.open_file(fa_name)

    for rec in SeqIO.parse(fah, "fasta"):
        seq_info[rec.id] = len(rec.seq)
        print rec.id, len(rec.seq)
    fah.close()
    
    print 
    print 'Number of FASTA entries: ', len(seq_info)
    for long_one in sorted(seq_info.items(), key=itemgetter(1), reverse=True):
        print 'Long contig length (bp): ', long_one[0], long_one[1]
        break
    for short_one in sorted(seq_info.items(), key=itemgetter(1)):
        print 'Short contig length (bp): ', short_one[0], short_one[1]
        break
    flength = 0 
    for ele in sorted(seq_info.items(), key=itemgetter(1)):
        flength += ele[1]
    print 'Average length of FASTA contig (bp): ', (flength/len(seq_info))
    print 


if __name__ == "__main__":
    print __doc__
    
