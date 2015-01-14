from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def U_to_T_nts(fasta_in, fasta_out):
    """
    change U to T 
    """

    fout = open(fasta_out, "w")
    for rec in SeqIO.parse(fasta_in, "fasta"):
    mut_seq = rec.seq.tomutable()

    for idx, nt in enumerate(mut_seq):
        if nt == 'U':
        mut_seq[idx]="T"

    mut_seq = SeqRecord(mut_seq, id=rec.id, description=rec.description)
    fasta_out.write(mut_seq.format("fasta"))
    fasta_out.close()

