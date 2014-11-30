#!/usr/bin/env python
"""
Fetch the details about the features explained in a GFF type file.

Usage: python feature_info.py in.gff

Requirements:
    gfftools : https://github.com/vipints/genomeutils/blob/master/gfftools

    TODO 
    read alignment max intron length 
    transcript prediction max exon length 
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

            if not ent1['exons'][idx].any():
                continue

            exon_cnt = len(ent1['exons'][idx])
            if exon_cnt > 1:

                intron_start = 0 
                for xq, excod in enumerate(ent1['exons'][idx]): 
                    
                    if xq > 0: 
                        #print intron_start, excod[0]-1 
                        if excod[0]-intron_start==1:
                            continue
                        # intron size 
                        intron_size[excod[0]-intron_start] = 1 
                        #print tid, excod[0]-intron_start

                    intron_start = excod[1]+1
                    exon_size[intron_start-excod[0]] = 1
    
    # sort the intron_size based on the keys 
    if intron_size:
        keys_int = sorted(intron_size)
        print 'MinIntronLength', int(keys_int[0]), int(keys_int[1]), int(keys_int[2])
        print 'MaxIntronLength', int(keys_int[-1]), int(keys_int[-2]), int(keys_int[-3])
        print 
        keys_ex = sorted(exon_size)
        print 'MinExonLength', int(keys_ex[0]), int(keys_ex[1]), int(keys_ex[2]) 
        print 'MaxExonLength', int(keys_ex[-1]), int(keys_ex[-2]), int(keys_ex[-3]) 
    else:
        print "Error in feature mapping, please check the source of parent child features" 
        print "May be the sources are different for parents and child features of the parent Gene"

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
