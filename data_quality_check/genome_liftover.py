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

    chain_file = ""

    # wrapping the liftover tool 
    genome_v_from = "hg16"
    genome_v_to = "hg19"

    #import ipdb 
    #ipdb.set_trace()


def data_processing(training_example_file):
    """
    read the ARTS supplement material to extract the example coordinates 
    """
    
    out_file = "hg16_examples.bed"
    out_file_fh = open(out_file, "w") 

    ## load arts v1 training dataset 
    for rec in SeqIO.parse(training_example_file, "fasta"):

        loc = re.match(r'(\w{4,5})([+-]?)(.*)', rec.id) ## chr10+1245980
        out_file_fh.write("%s\t%s\t%d\n" % (loc.group(1), loc.group(3), int(loc.group(3))+1))

    out_file_fh.close()
    return out_file


if __name__=="__main__":
    
    try:
        arts_training_data = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    chrom_bed_file = data_processing(arts_training_data) 

    lift_genome(chrom_bed_file)
