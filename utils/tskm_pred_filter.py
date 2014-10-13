#!/usr/bin/env python 
"""
Filter genes predicted by TranscriptSkimmer program  

based on the predicted splice-site consensus into account and length of the ORF 

Usage: python tskm_pred_filter.py in.gff in.fasta > filter.gff 
"""

import sys 
import numpy 
import collections
from Bio import SeqIO 
from gfftools import GFFParser

def __main__():

    try:
        gff_name = sys.argv[1]
        fas_file = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1) 

    gff_content = GFFParser.Parse(gff_name)
    
    # getting the spliced transcripts from the predicted gene list 
    transcripts_region = collections.defaultdict(list)
    for gene_recd in gff_content:
        spliced_transcript = collections.defaultdict(list)

        for idx, sub_rec in enumerate(gene_recd['transcripts']):
            exon_cnt = len(gene_recd['exons'][idx])

            if exon_cnt > 1: 
                for idk, ex in enumerate(gene_recd['exons'][idx]):
                    if idk == 0:
                        ex[0] = None 
                    if exon_cnt-1 == idk:
                        ex[1] = None

                    spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append(ex)

        transcripts_region[gene_recd['chr']].append(spliced_transcript)

    print len(transcripts_region)
    
    # check for consensus sequence 
    for fas_rec in SeqIO.parse(fas_file, "fasta"):
        if fas_rec.id in transcripts_region:
            for details in transcripts_region[fas_rec.id]:
                for genes, regions in details.items():

                    print genes 
                        
                    for region in regions:
                        print region

                        if genes[-1] == '+':
                            if not numpy.isnan(region[0]):
                                acc_seq = fas_rec.seq[int(region[0])-3:int(region[0])-1]
                                print acc_seq
                            if not numpy.isnan(region[1]):
                                don_seq = fas_rec.seq[int(region[1]):int(region[1])+2]
                                print don_seq  

                        elif genes[-1] == '-':
                            if not numpy.isnan(region[0]):
                                don_seq = fas_rec.seq[int(region[0])-3:int(region[0])-1]
                                don_seq = don_seq.reverse_complement()
                                print don_seq
                            
                            if not numpy.isnan(region[1]):
                                acc_seq = fas_rec.seq[int(region[1]):int(region[1])+2]
                                acc_seq = acc_seq.reverse_complement()
                                print acc_seq







    """
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
    """

if __name__ == "__main__":
    __main__() 
