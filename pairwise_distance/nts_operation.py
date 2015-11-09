#!/usr/bin/env python 
"""
function to change the nucleotide from a FASTA file
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def U_to_T_nts(fasta_in, fasta_out):
    """
    change nucleotide U to T 

    @args fasta_in: fasta file with "U" as a nucleotide 
    @type fasta_in: str
    @args fasta_out: result fasta file
    @type fasta_out: str 
    """

    fout = open(fasta_out, "w")

    for rec in SeqIO.parse(fasta_in, "fasta"):
        mut_seq = rec.seq.tomutable()

        for idx, nt in enumerate(mut_seq):
            if nt.upper() == 'U':
                mut_seq[idx]="T"
            else:
                mut_seq[idx] = nt.upper()

        mut_seq = SeqRecord(mut_seq, id=rec.id, description=rec.description)
        fasta_out.write(mut_seq.format("fasta"))

    fasta_out.close()
