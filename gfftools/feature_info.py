#!/usr/bin/env python
"""
Fetch the details about the features explained in a GFF type file.

example:
    MaxIntronLength MinIntronLength

Usage: python feature_info.py in.gff
"""

import re
import sys
import GFFParser

def Intron_det(TDB):
    """
    get the intron feature details MaxIntronLength MinIntronLength
    """

    intron_size = dict() 
    exon_size = dict() 
    for ent1 in TDB:
        for idx, tid in enumerate(ent1['transcripts']):

            exon_cnt = len(ent1['exons'][idx])
            if exon_cnt > 1:

                intron_start = 0 
                for xq, excod in enumerate(ent1['exons'][idx]): 
                    
                    if xq > 0: 
                        #print intron_start, excod[0]-1 
                        # intron size 
                        intron_size[excod[0]-intron_start] = 1 

                    intron_start = excod[1]+1
                    exon_size[intron_start-excod[0]] = 1
    
    # sort the intron_size based on the keys 
    keys_int = sorted(intron_size)
    print 'MinIntronLength', int(keys_int[0])
    print 'MaxIntronLength', int(keys_int[-1])
    print 
    keys_ex = sorted(exon_size)
    print 'MinExonLength', int(keys_ex[0])
    print 'MaxExonLength', int(keys_ex[-1])

def __main__():

    try:
        query_file = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    # get the annotated transcripts 
    Transdb = GFFParser.Parse(query_file)  

    # extract different features 
    Intron_det(Transdb)

if __name__ == "__main__": 
    __main__() 
