#!/usr/bin/env python 
"""
Make gff exon line for missing exon features of a gene defined in a GFF file using UTR's and CDS region

Usage: python make_gff_exon_line.py in.gff > out_exon.gff3  

Requirement:
    GFFParser : https://github.com/vipints/genomeutils/blob/master/gfftools/GFFParser.py

"""

import re
import sys
import GFFParser

def print_exon_line(tinfo):
    """
    writing exon feature info in GFF form 
    """

    for ent1 in tinfo:
        for idx, tid in enumerate(ent1['transcripts']): # - transcripts
            for idz, ex_cod in enumerate(ent1['exons'][idx]): # - exons 

                out_print = [ent1['chr'],
                            ent1['source'],
                            'exon',
                            str(int(ex_cod[0])),
                            str(int(ex_cod[1])),
                            '.',
                            ent1['strand'], 
                            '.',
                            'Parent=%s' % tid[0]]

                print '\t'.join(out_print) # - exon line  

            #break
        #break 


def __main__():

    try:
        query_file = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    transdb = GFFParser.Parse(query_file)  
    print_exon_line(transdb)


if __name__ == "__main__": 
    __main__() 
