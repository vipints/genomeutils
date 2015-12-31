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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict 



def lift_genome(chrom_bed_file):
    """
    run genome lifting program using a chain file 
    the chain file can be downloaded from ucsc main site
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


def write_lifted_example_seq(lifted_coord_file, genome_seq): 
    """
    TSS example sequence writer 
    """

    ## creating the features db 
    example_db = defaultdict(list) 

    fcod = open(lifted_coord_file, "rU")
    for line in fcod:
        line = line.strip('\n\r')
        parts = line.split("\t") 

        example_db[parts[0]].append((parts[0], parts[3], parts[1], parts[4], parts[5]))
    fcod.close() 
    ## ('chr1', '+', '155838683', '-1', 'NM_152280')

    out_pos_fh = open("tss_sig_pos_example.fa", 'w')
    out_neg_fh = open("tss_sig_neg_example.fa", 'w')

    boundary = 1200 
    ## getting the genome sequence 
    for rec in SeqIO.parse(genome_seq, "fasta"):
        if rec.id in example_db:
            for feat_det in example_db[rec.id]:

                motif_seq = rec.seq[int(feat_det[2])-boundary:int(feat_det[2])+boundary] # TSS ---A--

                if feat_det[1] == "-":
                    motif_seq = motif_seq.reverse_complement()
                
                # sanity check for the fetched sequence 
                if len(motif_seq) != boundary*2: 
                    continue
                if not motif_seq:
                    continue
                if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                    continue

                # write to fasta out 
                if feat_det[3] == "+1":

                    fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, feat_det[1], int(feat_det[2])), description='+1 %s' % feat_det[4])
                    out_pos_fh.write(fseq.format("fasta"))
                elif feat_det[3] == "-1":

                    fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, feat_det[1], int(feat_det[2])), description='-1 %s' % feat_det[4])
                    out_neg_fh.write(fseq.format("fasta"))

    
    out_pos_fh.close()
    out_neg_fh.close()



if __name__=="__main__":
    
    try:
        arts_training_data = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    chrom_bed_file = data_processing(arts_training_data) 
    lift_bed_file = lift_genome(chrom_bed_file)

    #lift_bed_file = "hg19_examples.bed"

    genome_file = "hg19.fa" 

    write_lifted_example_seq(lift_bed_file, genome_file)
    
