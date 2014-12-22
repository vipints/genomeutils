#!/usr/bin/python 
import sys, re
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord

try:
    fas_in=open(sys.argv[1], "rU")
except:
    sys.exit(-1)

seq_rec=[]
for rec in SeqIO.parse(fas_in, 'fasta'):
    seq_rec.append(str(rec.seq))
    seq_rec.append('N'*10000)
fas_in.close()
seq_rec.insert(0, 'N'*10000)
seq_rec = ''.join(seq_rec)

fasta_out = open("out.fa", 'w')
mut_seq = Seq(seq_rec)
seq_don = SeqRecord(mut_seq, id='ensemblGRCh37_biomart', description='hs-biomart-ensemblGRCh37')
fasta_out.write(seq_don.format("fasta"))
fasta_out.close()
