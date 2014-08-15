#!/usr/bin/env python 
"""
Filter genes predicted by TranscriptSkimmer program  

Usage: python tskm_pred_filter.py
"""

import sys 
import collections
from gfftools import GFFParser
from bx.intervals.cluster import ClusterTree


def __main__():

    try:
        gff_name = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1) 

    
    cluster_distance = 5 
    min_entries = 1 

    gff_content = GFFParser.Parse(gff_name)

    cluster_trees = collections.defaultdict(lambda:ClusterTree(cluster_distance, min_entries))

    feat_id_map = dict() 

    for xp, rec in enumerate(gff_content):
        feat_id_map[xp] = rec['name']

        cluster_trees[rec['chr']].insert(rec['start'], rec['stop'], xp)

    clean_entries = dict() 
    for chrom, sub_trees in cluster_trees.items():
        for start, stop, id in sub_trees.getregions():
            for xq in range(1):

                try:
                    clean_entries[feat_id_map[id[xq]]] = 0 
                except:
                    pass 
 
    for rec in gff_content:
        if rec['name'] in clean_entries:
            print '%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s' % (rec['chr'], rec['source'], rec['start'], rec['stop'], rec['strand'], rec['name'], rec['name'])
            for idx, tid in enumerate(rec['transcripts']):
                t_start = rec['exons'][idx][0][0]
                t_stop = rec['exons'][idx][-1][-1]
                print '%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s' % (rec['chr'], rec['source'], rec['transcript_type'][idx], t_start, t_stop, rec['strand'], tid[0], rec['name'])
                
                for ex_cod in rec['utr5_exons'][idx]:
                    print '%s\t%s\tfive_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s' % (rec['chr'], rec['source'], ex_cod[0], ex_cod[1], rec['strand'], tid[0])

                for ex_cod in rec['cds_exons'][idx]:
                    print '%s\t%s\tCDS\t%d\t%d\t.\t%s\t%d\tParent=%s' % (rec['chr'], rec['source'], ex_cod[0], ex_cod[1], rec['strand'], ex_cod[2], tid[0])

                for ex_cod in rec['utr3_exons'][idx]:
                    print '%s\t%s\tthree_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s' % (rec['chr'], rec['source'], ex_cod[0], ex_cod[1], rec['strand'], tid[0])

                for ex_cod in rec['exons'][idx]:
                    print '%s\t%s\texon\t%d\t%d\t.\t%s\t.\tParent=%s' % (rec['chr'], rec['source'], ex_cod[0], ex_cod[1], rec['strand'], tid[0])


if __name__ == "__main__":
    __main__() 
