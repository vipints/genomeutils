#!/usr/bin/env python 
"""
Filter genes predicted by TranscriptSkimmer program  

Usage: python tskm_pred_filter.py in.gff in.bam > filter.gff 
"""

import sys 
import pysam 
import collections
from gfftools import GFFParser
from bx.intervals.cluster import ClusterTree


def __main__():

    try:
        gff_name = sys.argv[1]
        bam_file = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1) 

    cluster_distance = 5 
    min_entries = 1 

    gff_content = GFFParser.Parse(gff_name)

    cluster_trees = collections.defaultdict(lambda:ClusterTree(cluster_distance, min_entries))

    feat_id_map = dict() 
    for xp, rec in enumerate(gff_content):

        ## bad fix to include the strand information 
        xp += 1 
        if rec['strand'] == '+':
            xp = str(xp) + '%d' % 0 
        elif rec['strand'] == '-':
            xp = str(xp) + '%d' % 1 

        feat_id_map[xp] = rec['name']
        cluster_trees[rec['chr']].insert(rec['start'], rec['stop'], int(xp))

    bam_fh = pysam.Samfile(bam_file, "rb")

    clean_entries = dict() 
    for chrom, sub_trees in cluster_trees.items():
        for start, stop, id in sub_trees.getregions():
            if len(id) > 1:

                reverse_cnt = 0 
                forward_cnt = 0 
                for read in bam_fh.fetch(str(chrom), start, stop):
                    if read.is_proper_pair and read.is_read1:
                        if read.is_reverse:
                            reverse_cnt +=1 
        
                        else:
                            forward_cnt +=1 
                
                ## Decision which transcript want to take
                ## --------------------------------------> gene 
                ## --->  <--- = 16    --->   <--- = 3000  = read count 
                ##  1      2           2       1         first and second read 

                selected = '1' # minus 
                selected = '0' if reverse_cnt > forward_cnt else  selected

                for element in id:
                    element=str(element)
                    xq = element[-1]
                    if xq == selected:
                        clean_entries[feat_id_map[element]] = 0 

            else: 
                clean_entries[feat_id_map[str(id[0])]] = 0 
    bam_fh.close()

    for rec in gff_content:
        if rec['name'] in clean_entries:
            print '%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s' % (rec['chr'], 
                    rec['source'], 
                    rec['start'], 
                    rec['stop'], 
                    rec['strand'], 
                    rec['name'], 
                    rec['name'])

            for idx, tid in enumerate(rec['transcripts']):
                t_start = rec['exons'][idx][0][0]
                t_stop = rec['exons'][idx][-1][-1]
                print '%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s' % (rec['chr'], 
                        rec['source'], 
                        rec['transcript_type'][idx], 
                        t_start, t_stop, 
                        rec['strand'], 
                        tid[0], 
                        rec['name'])
                
                for ex_cod in rec['utr5_exons'][idx]:
                    print '%s\t%s\tfive_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s' % (rec['chr'], 
                        rec['source'], 
                        ex_cod[0], 
                        ex_cod[1], 
                        rec['strand'], 
                        tid[0])

                for ex_cod in rec['cds_exons'][idx]:
                    print '%s\t%s\tCDS\t%d\t%d\t.\t%s\t%d\tParent=%s' % (rec['chr'], 
                        rec['source'], 
                        ex_cod[0], 
                        ex_cod[1], 
                        rec['strand'], 
                        ex_cod[2], 
                        tid[0])

                for ex_cod in rec['utr3_exons'][idx]:
                    print '%s\t%s\tthree_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s' % (rec['chr'], 
                        rec['source'], 
                        ex_cod[0], 
                        ex_cod[1], 
                        rec['strand'], 
                        tid[0])

                for ex_cod in rec['exons'][idx]:
                    print '%s\t%s\texon\t%d\t%d\t.\t%s\t.\tParent=%s' % (rec['chr'], 
                        rec['source'], 
                        ex_cod[0], 
                        ex_cod[1], 
                        rec['strand'], 
                        tid[0])


if __name__ == "__main__":
    __main__() 
