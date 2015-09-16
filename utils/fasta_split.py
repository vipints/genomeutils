#!/usr/bin/env python 

from Bio import SeqIO 


def split_fasta_records(fas_fname):
    """
    split the fasta records to multiple fasta file 
    """
    
    cnt = 0 
    records = 100 
    out_file_prefix = "tss_sig_pos_example_"

    for rec in SeqIO.parse(fas_fname, "fasta"):
     
        if cnt % records == 0:
            fname = "%s%d.fa" % (out_file_prefix, cnt) 

            try:
                fh.close()
            except:
                pass 
            
            fh = open(fname, "w")
            fh.write(rec.format("fasta"))
        else:
            fh.write(rec.format("fasta"))
                    
        cnt +=1 


if __name__ == "__main__":
    print __doc__
    
