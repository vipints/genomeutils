#!/usr/bin/env
"""
"""
import sys, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    fasta_in=sys.argv[1]
except:
    print __doc__
    sys.exit(-1)

fh=open(fasta_in, 'rU')
for rec in SeqIO.parse(fh, 'fasta'):
    print rec.id 
    sep=500*'N'
    rec.seq=rec.seq.split(sep)

    for idx, record in enumerate(rec.seq):
        if not record:
            continue 
        fseq = SeqRecord(record, id='%s_%d' % (rec.id, idx), description='')

        out_min_fh = open("%s_%d.fa" % (rec.id, idx), "w") 
        out_min_fh.write(fseq.format("fasta"))
        out_min_fh.close()

        #break
    #print len(rec.seq)
    break
fh.close()


