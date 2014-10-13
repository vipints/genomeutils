#!/usrt/bin/env
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

def __main__():
    
    try:
        arts_labels = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    genome_v_from = "hg16"
    genome_v_to = "hg19"

    # fetch the genome coordinates 
    for rec in SeqIO.parse(arts_labels, "fasta"):
        id = re.split("\+|\-", rec.id) 
        
        print '%s:%s-%s' % (id[0], id[1], id[1])

        #break

    # wrapping the liftover tool 



if __name__=="__main__":
    __main__() 
