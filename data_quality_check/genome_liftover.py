#!/usr/bin/env python 
"""
fetch genome coordinates from old version of human genome and run UCSC liftover tool based on the genome version

    ARTS labels data with fasta sequence
    - liftover program  
    - genome versions 
    - fasta seq
"""

import re 
import sys
from Bio import SeqIO 



def lift_genome(chrom_bed_file):
    """
    run genome lifting program using a chain file 
    """

    chain_file = "hg16ToHg19.chain"
    # TODO wrapping the liftover tool 
    lifted_file = "hg19_examples.bed" 

    return lifted_file



def data_processing(training_example_file):
    """
    read the ARTS supplement material to extract the example coordinates 

    creating a bed file  
    """
    
    out_file = "hg16_examples.bed"
    out_file_fh = open(out_file, "w") 

    ## load arts v1 training dataset 
    for rec in SeqIO.parse(training_example_file, "fasta"):

        loc = re.match(r'(\w{4,5})([+-]?)(.*)', rec.id) ## chr10+1245980
        desc = rec.description.split(" ") ## chr9-178721 +1 NM_6587384

        out_file_fh.write("%s\t%s\t%d\t%s\t%s\t%s\n" % (loc.group(1), loc.group(3), int(loc.group(3))+1, loc.group(2), desc[1], desc[2]))

    out_file_fh.close()
    return out_file


def write_lifted_example_seq(): 
    """
    TSS example sequence writer 
    """
    
    lifted_coord_file = ""
    fcod = open(lifted_coord_file, "rU")
     
    for rec in fcod:
        rec = rec.strip('\n\r')

    fcod.close() 

    out_pos_fh = open("tss_sig_pos_example.fa", 'w')
    out_neg_fh = open("tss_sig_neg_example.fa", 'w')

    foh = helper.open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():

                    #motif_seq = rec.seq[int(loc[0])-boundary:int(loc[0])+boundary+1] # TSS ---A---
                    motif_seq = rec.seq[int(loc[0])-boundary:int(loc[0])+boundary] # TSS ---A--

                    if loc[1] == '-': 
                        motif_seq = motif_seq.reverse_complement()

                    # sanity check for the fetched sequence 
                    if len(motif_seq) != boundary*2: 
                        continue
                    if not motif_seq:
                        continue
                    if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                        continue

                    # write to fasta out 
                    fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, loc[1], int(loc[0])), description='+1 %s' % fid)
                    out_pos_fh.write(fseq.format("fasta"))

    out_pos_fh.close()
    out_neg_fh.close()

    foh.close()
    


if __name__=="__main__":
    
    """
    try:
        arts_training_data = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    """

    #chrom_bed_file = data_processing(arts_training_data) 
    #lift_bed_file = lift_genome(chrom_bed_file)

    lift_bed_file = "hg19_examples.bed"
    genome_file = "" 

    #write_lifted_example_seq()
    
