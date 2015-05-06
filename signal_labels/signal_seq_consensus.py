#!/usr/bin/env python
"""
Consensus sequence specificity at Translational Initiation Site [TIS], 
Translation Stop Site [cdsStop] and  splice sites [don, acc] signals 
from the genome annotation and genome sequence. 

Usage: python signal_seq_consensus.py in.fasta.(gz) in.gtf.(gz)

Requirements:
    BioPython:- http://biopython.org 
    gfftools:- https://github.com/vipints/gfftools
"""


from __future__ import division
import re
import os 
import sys
import random
import numpy as np 
from Bio import SeqIO 
from collections import defaultdict
from gfftools import helper, GFFParser 


def check_consensus(gfname, faname):
    """
    checking the consensus sequence of the specific signal 

    @args gfname: genome annotation in gtf/gff file
    @type gfname: str  
    @args faname: genome sequence in fasta file 
    @type faname: str  
    """

    ## extract genome annotation from gtf file 
    gtf_file_content = GFFParser.Parse(gfname)
    print 'processed annotation file'

    ## signals considering 
    for signal in ['splice', 'tis']:
        
        gtf_db, feature_cnt = get_label_regions(gtf_file_content, signal)
        print 'extracted %d %s signal regions' % (feature_cnt, signal) 

        if signal == 'splice':

            don_true_seq, acc_true_seq, don_fal_seq, acc_fal_seq = true_ss_seq_fetch(faname, gtf_db, boundary=100)
            print 'summary of', signal, 'signal consensus'

            print 'don site cons %f non-cons %f' % (round((don_true_seq/feature_cnt)*100, 2), round((don_fal_seq/feature_cnt)*100, 2))
            print 'acc site cons %f non-cons %f' % (round((acc_true_seq/feature_cnt)*100, 2), round((acc_fal_seq/feature_cnt)*100, 2))

            break ## one singal 


def get_label_regions(gtf_content, signal):
    """
    get signal sequence location from the annotation

    @args gtf_content: parsed genome annotation from GFF/GTF 
    @type gtf_content: numpy.array object 
    @args signal: genome signal type 
    @type signal: str 
    """

    feat_cnt = 0
    anno_db = defaultdict(list) 
    
    for feature in gtf_content: # gene list returned from GFFParse function 
        mod_anno_db = dict()

        if signal == 'tis':
            for xp, ftid in enumerate(feature['transcripts']):
                if feature['cds_exons'][xp].any():
                    feat_cnt += 1
                    mod_anno_db[ftid[0]] = (feature['tis'][xp],
                                feature['strand'],
                                (int(feature['start']), int(feature['stop']))
                                )
        elif signal == 'cdsstop':
            for xp, ftid in enumerate(feature['transcripts']):
                if feature['cds_exons'][xp].any():
                    feat_cnt += 1
                    mod_anno_db[ftid[0]] = (feature['cdsStop'][xp],
                                feature['strand'],
                                (int(feature['start']), int(feature['stop']))
                                )
        elif signal == 'splice':
            for xp, ftid in enumerate(feature['transcripts']):
                if np.isnan(feature['exons'][xp]).any():
                    continue #clearing empty entiies 

                if len(feature['exons'][xp]) > 2: 
                    id_cnt = 1
                    for ex in feature['exons'][xp][1:-1]:
                        feat_cnt += 1
                        alt_id = ftid[0] + '.' + str(id_cnt)
                        id_cnt += 1
                        mod_anno_db[alt_id] = (ex[0], 
                                            ex[1], 
                                            feature['strand'])
        if mod_anno_db:
            anno_db[feature['chr']].append(mod_anno_db)
    
    return dict(anno_db), feat_cnt 

def true_ss_seq_fetch(fnam, Label, boundary):
    """
    true splice signals 
    """

    foh = helper.open_file(fnam)

    don_cnt_pl = don_cnt_mi = acc_cnt_pl = acc_cnt_mi = 0 
    don_in_pl = don_in_mi = acc_in_pl = acc_in_mi = 0 

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lfeat in Label[rec.id]:
                for fid, loc in Lfeat.items():
                
                    acc_ind = don_ind = 0 
                    if loc[-1] == '+': 
                        acc_mot_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]
                        if not acc_mot_seq:
                            acc_ind = 1
                        if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                            acc_ind = 1
                   
                        if acc_ind:
                            acc_in_pl += 1
                        else:
                            acc_cnt_pl += 1 

                        don_mot_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                        if not don_mot_seq:
                            don_ind = 1 
                        if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                            don_ind = 1 

                        if don_ind:
                            don_in_pl += 1 
                        else:
                            don_cnt_pl += 1 

                    elif loc[-1] == '-':
                        don_mot_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]

                        if not don_mot_seq:
                            don_ind = 1 
                        if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'AC':
                            don_ind = 1 

                        if don_ind:
                            don_in_mi += 1 
                        else:
                            don_cnt_mi += 1 

                        acc_mot_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                        if not acc_mot_seq:
                            acc_ind = 1
                        if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'CT':
                            acc_ind = 1

                        if acc_ind:
                            acc_in_mi += 1 
                        else:
                            acc_cnt_mi += 1 

    print 
    print '%d and %d Consensus DONOR sites on positive and negative strand' % (don_cnt_pl, don_cnt_mi) 
    print '%d and %d Non-Consensus DONOR sites on  positive and negative strand' % (don_in_pl, don_in_mi)
    print 
    print '%d and %d Consensus ACCEPTOR sites on positive and negative strand' % (acc_cnt_pl, acc_cnt_mi) 
    print '%d and %d Non-Consensus ACCEPTOR sites on positive and negative strand' % (acc_in_pl, acc_in_mi)
    print 
    foh.close()
    return don_cnt_pl+don_cnt_mi, acc_cnt_pl+acc_cnt_mi, don_in_pl+don_in_mi, acc_in_pl+acc_in_mi

if __name__ == "__main__":

    try:
        faname = sys.argv[1] # sequence 
        gfname = sys.argv[2] # annotation 
    except:
        print __doc__
        sys.exit(-1)
    
    check_consensus(gfname, faname) 
