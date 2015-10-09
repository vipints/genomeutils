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



def lift_genome(trainig_example_file):
    """
    """
    
    ## load arts v1 training dataset 
    for rec in SeqIO.parse(trainig_example_file, "fasta"):

        location = re.match(r'(\w{1,2})([+-]?)(.*)', rec.id)

        

        break

    # wrapping the liftover tool 
    genome_v_from = "hg16"
    genome_v_to = "hg19"



if __name__=="__main__":
    
    try:
        arts_training_data = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    lift_genome(arts_training_data)
