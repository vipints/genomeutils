#!/usr/bin/env python
"""
Extract genome signal sequence labels like Transcription Start Site [TSS], 
Translational Initiation Site [TIS] and Splice Sites [don, acc] from the 
genome annotation in GTF/GFF with genome sequence in FASTA format. 

The number of positive and negative labels for each signal sequence are 
generated randomly from each chromosome. The numbers can be adjusted. 
Currently which is set to TT and NN for positive and negative respectively. 

The positive labels are generated from YY nucleotides upstream and downstream 
of the annotated genome signal region. Extracted labels are stored in the 
base folder of input file with signal specific names. 
ex: TSS [tss_sig_{minus|plus}_label.fa]

Usage: python generate_genome_seq_labels.py in.fasta.(gz|bz2) in.gtf.(gz|bz2)

Requirements:
    BioPython:- http://biopython.org 
    gfftools:- https://github.com/vipints/genomeutils/tree/master/gfftools
"""

import re
import os 
import sys
import random
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from gfftools import helper, GFFParser 

def __main__():
    """
    Main processing unit
    """
    try:
        faname = sys.argv[1] # sequence 
        gfname = sys.argv[2] # annotation 
    except:
        print __doc__
        sys.exit(-1)

    # adjust the training label sequence count & flaking region length  
    label_cnt = 2 # number of labels 
    sp_boundary = 100 # flanking region nucleotides to the splice site 
    t_boundary = 200 # flanking region nucleotides to TIS and TSS 

    # extract genome annotation from gtf/gff type file 
    anno_file_content = GFFParser.Parse(gfname)
    print 'processed annotation file'
    
    # genomic signals : don/acc - Transcription - Translation 
    for signal in ['splice', 'tss', 'tis']: # 
        
        gtf_db, feature_cnt = get_label_regions(anno_file_content, signal)
        print 'extracted', feature_cnt, signal, 'signal regions'

        posLabel, COUNT = select_labels(gtf_db, feature_cnt, label_cnt) 
        print 'selecting', COUNT, 'random', signal, 'labels'

        if signal == 'splice':
            true_ss_seq_fetch(faname, posLabel, sp_boundary) 
            print 'fetched don/acc plus signal lables'

            false_ss_seq_fetch(faname, posLabel, sp_boundary)
            print 'fetched don/acc minus signal lables'

        elif signal == 'tis':
            true_tis_seq_fetch(faname, posLabel, t_boundary)
            print 'fetched', signal, 'plus signal lables'

            false_tis_seq_fetch(faname, posLabel, t_boundary)
            print 'fetched', signal, 'minus signal lables'

        else:
            plus_tss_seq_fetch(faname, posLabel, t_boundary)
            print 'fetched', signal, 'plus signal lables'

            minus_tss_seq_fetch(faname, posLabel, t_boundary)
            print 'fetched', signal, 'minus signal lables'
        
        #TODO remove the extra features 

def false_tis_seq_fetch(fnam, Label, boundary):
    """
    fetch the minus TIS signal sequence
    """
    foh = helper._open_file(fnam)
    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    #out_pos_fh = open(out_path + "/" + "tis_sig_plus_label.fa", 'w')
    out_min_fh = open("tis_sig_minus_label.fa", 'w')

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    if loc[-1] == '+': 
                        motif_seq = rec.seq[int(loc[0])-boundary:int(loc[0])+boundary+1]
                        # get index for negative signal label sequence site 
                        idx = [xq.start() for xq in re.finditer(re.escape('ATG'), str(motif_seq).upper())]
                        # removing the true signl sequence site from selected false sites
                        if boundary-1 in idx:
                            idx.remove(boundary-1)
                        if not idx:
                            continue 
                        # limit to take maximum 2 false labels from one defined feature
                        if len(idx) > 2:
                            idx = random.sample(idx, 2)
                        # get the false labels for randomly selected region or the available ones   
                        for xp in idx:
                            # adjusting the coordinate to the false site 
                            rloc_min = (int(loc[0])-boundary)+xp
                            motif_seq = rec.seq[(rloc_min-boundary)+1:(rloc_min+boundary)+2]
                            # check for sanity and consensus of the fetched sequence region 
                            if not motif_seq:
                                continue
                            if 'N' in motif_seq.upper():
                                continue
                            if str(motif_seq[boundary-1:boundary+2]).upper() != 'ATG':
                                continue
                            # result to fasta out
                            fseq = SeqRecord(motif_seq.upper(), id=fid, description='-ve label')
                            out_min_fh.write(fseq.format("fasta"))

                    elif loc[-1] == '-': 
                        motif_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-1]
                        # get index for negative signal label sequence site 
                        idx = [xq.start() for xq in re.finditer(re.escape('CAT'), str(motif_seq).upper())]
                        # removing the true signal sequence site from selected false sites
                        if boundary-1 in idx:
                            idx.remove(boundary-1)
                        if not idx:
                            continue 
                        # limit to take maximum 2 false labels from one defined feature
                        if len(idx) > 2:
                            idx = random.sample(idx, 2)
                        # get the false labels for randomly selected region or the available ones   
                        for xp in idx:
                            # adjusting the coordinate to the false site 
                            rloc_min = (int(loc[0])-boundary)+xp
                            motif_seq = rec.seq[(rloc_min-boundary)-1:(rloc_min+boundary)]
                            motif_seq = motif_seq.reverse_complement()
                            # check for sanity and consensus of the fetched sequence region 
                            if not motif_seq:
                                continue
                            if 'N' in motif_seq.upper():
                                continue
                            if str(motif_seq[boundary-1:boundary+2]).upper() != 'ATG':
                                continue
                            # result to fasta out
                            fseq = SeqRecord(motif_seq.upper(), id=fid, description='-ve label')
                            out_min_fh.write(fseq.format("fasta"))
    out_min_fh.close()
    foh.close()

def true_tis_seq_fetch(fnam, Label, boundary):
    """
    fetch the plus TIS signal sequence 
    """
    foh = helper._open_file(fnam)
    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    #out_pos_fh = open(out_path + "/" + "tis_sig_plus_label.fa", 'w')
    out_pos_fh = open("tis_sig_plus_label.fa", 'w')

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    if loc[-1] == '+': 
                        motif_seq = rec.seq[int(loc[0])-boundary:int(loc[0])+boundary+1]
                        # check for sanity and consensus of the fetched sequence region 
                        if not motif_seq:
                            continue
                        if 'N' in motif_seq.upper():
                            continue
                        if str(motif_seq[boundary-1:boundary+2]).upper() != 'ATG':
                            continue
                        # result to fasta out
                        fseq = SeqRecord(motif_seq.upper(), id=fid, description='+ve label')
                        out_pos_fh.write(fseq.format("fasta"))

                    elif loc[-1] == '-': 
                        motif_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-1]
                        motif_seq = motif_seq.reverse_complement()
                        # check for sanity and consensus of the fetched sequence region 
                        #if len(motif_seq) != 2*boundary: 
                        #    continue
                        if not motif_seq:
                            continue
                        if 'N' in motif_seq.upper():
                            continue
                        if str(motif_seq[boundary-1:boundary+2]).upper() != 'ATG':
                            continue
                        # result to fasta out
                        fseq = SeqRecord(motif_seq.upper(), id=fid, description='+ve label')
                        out_pos_fh.write(fseq.format("fasta"))
    out_pos_fh.close()
    foh.close()
        
def get_label_regions(gtf_content, signal):
    """
    get signal sequence location from the annotation
    """
    feat_cnt = 0
    anno_db = defaultdict(list) 
    
    for feature in gtf_content: # gene list returned from GFFParse function 
        mod_anno_db = dict()
        if signal == 'tss':
            for xp, ftid in enumerate(feature['transcripts']):
                feat_cnt += 1
                mod_anno_db[ftid[0]] = (feature['tss'][xp], 
                                feature['strand'])
        elif signal == 'tis':
            for xp, ftid in enumerate(feature['transcripts']):
                if feature['cds_exons'][xp].any():
                    feat_cnt += 1
                    mod_anno_db[ftid[0]] = (feature['tis'][xp],
                                feature['strand'])
        elif signal == 'splice':
            # going through each transcripts annotated 
            for xp, ftid in enumerate(feature['transcripts']):
                # spliced transcripts 
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

def false_ss_seq_fetch(fnam, Label, boundary):
    """
    false splice signals
    """
    foh = helper._open_file(fnam)
    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam) 
    #don_min_fh = open(out_path + "/don_sig_minus_label.fa", 'w')
    #acc_min_fh = open(out_path + "/acc_sig_minus_label.fa", 'w')
    don_min_fh = open("don_sig_minus_label.fa", 'w')
    acc_min_fh = open("acc_sig_minus_label.fa", 'w')

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for fLabel in Label[rec.id]:
                for fid, loc in fLabel.items():
                    if loc[-1]=='+':
                        # acceptor splice site signal 
                        acc_t_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]
                        # get index for negative signal label sequence site 
                        idx = [xq.start() for xq in re.finditer(re.escape('AG'), str(acc_t_seq).upper())]
                        # removing the true signal sequence site from selected false sites
                        if boundary-1 in idx:
                            idx.remove(boundary-1)
                        if not idx:
                            continue 
                        # limit to take maximum 2 false labels from one defined feature
                        if len(idx) > 2:
                            idx = random.sample(idx, 2)
                        for xp in idx:
                            # adjusting the coordinate to the false site 
                            rloc_min = (int(loc[0])-boundary)+xp
                            acc_mot_seq = rec.seq[(rloc_min-boundary)-1:(rloc_min+boundary)-1]
                            # check for sanity and consensus of the fetched sequence region 
                            #if len(acc_mot_seq) != 2*boundary:
                            #    continue
                            if not acc_mot_seq:
                                continue
                            if 'N' in acc_mot_seq.upper():
                                continue
                            if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                                continue 
                            # write to fasta out 
                            fseq_acc = SeqRecord(acc_mot_seq.upper(), id=fid, description='-ve label')
                            acc_min_fh.write(fseq_acc.format("fasta"))
                        
                        # donor splice site signal 
                        don_t_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                        # get index for negative signal label sequence site 
                        idx = [xq.start() for xq in re.finditer(re.escape('GT'), str(don_t_seq).upper())]
                        # removing the true signal sequence site from selected false sites
                        if boundary-1 in idx:
                            idx.remove(boundary-1)
                        if not idx:
                            continue
                        # limit to take maximum 2 false labels from one defined feature
                        if len(idx) > 2:
                            idx = random.sample(idx, 2)
                        for xj in idx:
                            # adjusting the coordinate to the false site 
                            rloc_pos = (int(loc[1])-boundary)+xj
                            don_mot_seq = rec.seq[(rloc_pos-boundary)+2:(rloc_pos+boundary)+2]
                            # check for sanity and consensus of the fetched sequence region 
                            #if len(don_mot_seq) != 2*boundary:
                            #    continue
                            if not don_mot_seq:
                                continue
                            if 'N' in don_mot_seq.upper():
                                continue
                            if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                                continue 
                            # write to fasta out 
                            fseq_don = SeqRecord(don_mot_seq.upper(), id=fid, description='-ve label')
                            don_min_fh.write(fseq_don.format("fasta"))

                    elif loc[-1]=="-":
                        # donor splice signal site 
                        don_t_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]
                        # get index for negative signal label sequence site 
                        idx = [xq.start() for xq in re.finditer(re.escape('AC'), str(don_t_seq).upper())]
                        # removing the true signal sequence site from selected false sites
                        if boundary-1 in idx:
                            idx.remove(boundary-1)
                        if not idx:
                            continue
                        # limit to take maximum 2 false labels from one defined feature
                        if len(idx) > 2:
                            idx = random.sample(idx, 2)
                        for xp in idx:
                            # adjusting the coordinate to the false site 
                            rloc_min = (int(loc[0])-boundary)+xp
                            don_mot_seq = rec.seq[(rloc_min-boundary)-1:(rloc_min+boundary)-1]
                            don_mot_seq = don_mot_seq.reverse_complement()
                            # check for sanity and consensus of the fetched sequence region 
                            #if len(don_mot_seq) != 2*boundary:
                            #    continue 
                            if not don_mot_seq:
                                continue
                            if 'N' in don_mot_seq.upper():
                                continue
                            if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                                continue 
                            # write to fasta out 
                            fseq_don = SeqRecord(don_mot_seq.upper(), id=fid, description='-ve label')
                            don_min_fh.write(fseq_don.format("fasta"))
                        
                        # acceptor splice signal site 
                        acc_t_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                        # get index for negative signal label sequence site 
                        idx = [xq.start() for xq in re.finditer(re.escape('CT'), str(acc_t_seq).upper())]
                        # removing the true signal sequence site from selected false sites
                        if boundary-1 in idx:
                            idx.remove(boundary-1)
                        if not idx:
                            continue
                        # limit to take maximum 2 false labels from one defined feature
                        if len(idx) > 2:
                            idx = random.sample(idx, 2)
                        for xk in idx:
                            # adjusting the coordinate to the false site 
                            rloc_pos = (int(loc[1])-boundary)+xk
                            acc_mot_seq = rec.seq[(rloc_pos-boundary)+2:(rloc_pos+boundary)+2]
                            acc_mot_seq = acc_mot_seq.reverse_complement()
                            # check for sanity and consensus of the fetched sequence region 
                            #if len(acc_mot_seq) != 2*boundary:
                            #    continue
                            if not acc_mot_seq:
                                continue
                            if 'N' in acc_mot_seq.upper():
                                continue
                            if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                                continue 
                            # write to fasta out 
                            fseq_acc = SeqRecord(acc_mot_seq.upper(), id=fid, description='-ve label')
                            acc_min_fh.write(fseq_acc.format("fasta"))
    foh.close()
    don_min_fh.close()
    acc_min_fh.close()

def minus_tss_seq_fetch(fnam, Label, boundary):
    """
    fetch the minus TSS signal sequence label
    """
    foh = helper._open_file(fnam)
    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam) ## result to fasta base dir 
    #out_min_fh = open(out_path + "/" + "_sig_minus_label.fa", 'w')
    out_min_fh = open("tss_sig_minus_label.fa", 'w')

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    for ndr in range(3):
                        # max 3 labels from a true feature 
                        rloc = random.randint(int(loc[0])-boundary,int(loc[0])+boundary)
                        # remove the true signal index from random sampling 
                        if rloc == int(loc[0]):
                            continue

                        motif_seq = rec.seq[rloc-boundary:rloc+boundary+1]
                        # sanity check for featched sequence 
                        #if len(motif_seq) != 2*boundary:
                        #    continue
                        if not motif_seq:
                            continue
                        if 'N' in motif_seq.upper():
                            continue
                        # write to fasta out 
                        fseq = SeqRecord(motif_seq.upper(), id=fid+'_'+str(ndr), description='-ve label')
                        out_min_fh.write(fseq.format("fasta"))
    out_min_fh.close()
    foh.close()

def true_ss_seq_fetch(fnam, Label, boundary):
    """
    True splice signals 
    """
    foh = helper._open_file(fnam)
    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    #don_pos_fh = open(out_path + "/don_sig_plus_label.fa", 'w')
    #acc_pos_fh = open(out_path + "/acc_sig_plus_label.fa", 'w')
    don_pos_fh = open("don_sig_plus_label.fa", 'w')
    acc_pos_fh = open("acc_sig_plus_label.fa", 'w')

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lfeat in Label[rec.id]:
                for fid, loc in Lfeat.items():
                    if loc[-1] == '+': 
                        # acceptor splice site 
                        acc_mot_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]
                        # sanity check and consensus sequence from fetched sequence 
                        #if len(acc_mot_seq) != 2*boundary:
                        #    continue
                        if not acc_mot_seq:
                            continue
                        if 'N' in acc_mot_seq.upper():
                            continue
                        if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                            continue 
                        # write to fasta out 
                        fseq_acc = SeqRecord(acc_mot_seq.upper(), id=fid, description='+ve label')
                        acc_pos_fh.write(fseq_acc.format("fasta"))
                        # donor splice site 
                        don_mot_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                        #if len(don_mot_seq) != 2*boundary:
                        #    continue 
                        if not don_mot_seq:
                            continue
                        if 'N' in don_mot_seq.upper():
                            continue
                        if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                            continue 
                        # write to fasta out 
                        fseq_don = SeqRecord(don_mot_seq.upper(), id=fid, description='+ve label')
                        don_pos_fh.write(fseq_don.format("fasta"))

                    elif loc[-1] == '-':
                        # donor splice site signal 
                        don_mot_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]
                        don_mot_seq = don_mot_seq.reverse_complement()
                        # sanity check and consensus sequence from fetched sequence 
                        #if len(don_mot_seq) != 2*boundary:
                        #    continue 
                        if not don_mot_seq:
                            continue
                        if 'N' in don_mot_seq.upper():
                            continue
                        if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                            continue 
                        # write to fasta out 
                        fseq_don = SeqRecord(don_mot_seq.upper(), id=fid, description='+ve label')
                        don_pos_fh.write(fseq_don.format("fasta"))
                        # acceptor splice signal 
                        acc_mot_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                        acc_mot_seq = acc_mot_seq.reverse_complement()
                        # sanity check and consensus sequence from fetched sequence 
                        #if len(acc_mot_seq) != 2*boundary:
                        #    continue
                        if not acc_mot_seq:
                            continue
                        if 'N' in acc_mot_seq.upper():
                            continue
                        if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                            continue 
                        # write to fasta out 
                        fseq_acc = SeqRecord(acc_mot_seq.upper(), id=fid, description='+ve label')
                        acc_pos_fh.write(fseq_acc.format("fasta"))
    don_pos_fh.close()
    acc_pos_fh.close()
    foh.close()

def plus_tss_seq_fetch(fnam, Label, boundary):
    """
    fetch the plus TSS signal sequence 
    """
    foh = helper._open_file(fnam)
    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    #out_pos_fh = open(out_path + "/" + "_sig_plus_label.fa", 'w')
    out_pos_fh = open("tss_sig_plus_label.fa", 'w')

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    motif_seq = rec.seq[int(loc[0])-boundary:int(loc[0])+boundary+1]
                    if loc[-1] == '-': 
                        motif_seq = motif_seq.reverse_complement()
                    # sanity check for the fetched sequence 
                    #if len(motif_seq) != 2*boundary: 
                    #    continue
                    if not motif_seq:
                        continue
                    if 'N' in motif_seq.upper():
                        continue
                    # write to fasta out 
                    fseq = SeqRecord(motif_seq.upper(), id=fid, description='+ve label')
                    out_pos_fh.write(fseq.format("fasta"))
    out_pos_fh.close()
    foh.close()

def select_labels(feat_db, feat_count, label_cnt):
    """
    Random sampling to select signal lables
    """
    assert label_cnt <= feat_count, 'Number of features annotated ' + str(feat_count)

    try:
        accept_prob = (1.0*label_cnt)/feat_count
    except:
        accept_prob = 1

    while 1: # ensure the label count 
        counter, LSet = recursive_fn(feat_db, label_cnt, accept_prob)
        if label_cnt <= counter:
            break
        print '    recursive ...', counter

    return LSet, counter

def recursive_fn(f_db, lb_cnt, apt_prob):
    """
    This function returns the random samples based on the label counts 
    """
    pLabel = defaultdict(list)
    cnt = 0 
    for chrom, feat in f_db.items():
        tmp_db = dict()
        for sub_feat in feat:
            for fid, location in sub_feat.items():
                rnb = random.random()
                if lb_cnt == cnt:
                    break
                if rnb <= apt_prob:
                    cnt += 1
                    tmp_db[fid] = location
        pLabel[chrom].append(tmp_db)

    return cnt, dict(pLabel)

if __name__=="__main__":
    __main__()
