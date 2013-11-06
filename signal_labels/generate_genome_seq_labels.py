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

Usage: python generate_genome_seq_labels.py in.fasta.(gz) in.gtf.(gz)

Requirements:
    BioPython:- http://biopython.org 
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

    # look for a fasta index file in the base folder 
    #for fa_index_file in os.listdir(os.path.dirname(os.path.realpath(faname))):
    #    if fa_index_file.endswith(".fai"):
    #        break
    #print fa_index_file

    #contigs = dict()
    #faih = helper._open_file(os.path.dirname(os.path.realpath(faname)) + '/' + fa_index_file)
    #for chr in faih:
    #    chr = chr.strip().split('\t')
    #    contigs[chr[0]] = 1
    #print 'selected super contigs'

    # extract genome annotation from gtf/gff file 
    anno_file_content = GFFParser.Parse(gfname)
    print 'processed annotation file'
    
    # genomic signals
    for signal in ['splice', 'TIS', 'TSS']:

        gtf_db, feature_cnt = get_label_regions(anno_file_content, signal)
        print 'extracted', feature_cnt, signal, 'signal regions'

        posLabel, COUNT = select_labels(gtf_db, feature_cnt, label_cnt=4) # number of labels 
        print 'selecting', COUNT, 'random', signal, 'labels'

        if signal == 'splice':
            
            true_ss_seq_fetch(faname, posLabel, boundary=100) # flanking nucleotides 
            print 'fetched don/acc plus signal lables'
            false_ss_seq_fetch(faname, posLabel, boundary=100)
            print 'fetched don/acc minus signal lables'

        break 

    #    # TODO check the consensus of regions extracted from the following code

    #    else:
    #        #TODO one set random labels are fine for both plus and minus labels 
    #        pos_seq_fetch(faname, posLabel, signal, boundary=100)
    #        print 'fetched ' ,signal, ' plus signal lables'
    #        
    #        min_seq_fetch(faname, posLabel, signal, boundary=100)
    #        print 'fetched ' ,signal, ' minus signal lables'
        

def get_label_regions(gtf_content, signal):
    """
    get signal sequence location from the annotation
    """
    feat_cnt = 0
    anno_db = defaultdict(list) 
    
    for feature in gtf_content: # gene list in numpy format 
        #if not feature['chr'] in chrom:
        #    continue

        mod_anno_db = dict()
        if signal == 'TSS':
            for xp, ftid in enumerate(feature['transcripts']):
                feat_cnt += 1
                mod_anno_db[ftid[0]] = (feature['tss'][xp], 
                                feature['strand'])
        elif signal == 'TIS':
            for xp, ftid in enumerate(feature['transcripts']):
                if feature['cds_exons'][xp].any():
                    feat_cnt += 1
                    mod_anno_db[ftid[0]] = (feature['tis'][xp],
                                feature['strand'])
        elif signal == 'splice':
            for xp, ftid in enumerate(feature['transcripts']): # considering each transcript of a gene
                if len(feature['exons'][xp]) > 2: # spliced transcript 
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
    don_min_fh = open(out_path + "/don_sig_minus_label.fa", 'w')
    acc_min_fh = open(out_path + "/acc_sig_minus_label.fa", 'w')

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for fLabel in Label[rec.id]:
                for fid, loc in fLabel.items():
                    if loc[-1]=='+':
                        acc_t_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]
                        # index for negative sites 
                        idx = [xq.start() for xq in re.finditer(re.escape('AG'), str(acc_t_seq).upper())]
                        if boundary-1 in idx:
                            idx.remove(boundary-1)
                        if not idx:
                            continue 

                        for xp in idx:
                            rloc_min = (int(loc[0])-boundary)+xp
                            acc_mot_seq = rec.seq[(rloc_min-boundary)-1:(rloc_min+boundary)-1]
                            if len(acc_mot_seq) != 2*boundary:
                                continue
                            if not acc_mot_seq:
                                continue
                            if 'N' in acc_mot_seq.upper():
                                continue
                            if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                                continue 

                            fseq_acc = SeqRecord(acc_mot_seq.upper(), id=fid, description='-ve label')
                            acc_min_fh.write(fseq_acc.format("fasta"))
                            break
                        ##
                        don_t_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                        # indexing the negative sites 
                        idx = [xq.start() for xq in re.finditer(re.escape('GT'), str(don_t_seq).upper())]
                        if boundary-1 in idx:
                            idx.remove(boundary-1)
                        if not idx:
                            continue

                        for xj in idx:
                            rloc_pos = (int(loc[1])-boundary)+xj
                            don_mot_seq = rec.seq[(rloc_pos-boundary)+2:(rloc_pos+boundary)+2]
                    
                            if len(don_mot_seq) != 2*boundary:
                                continue 
                            if not don_mot_seq:
                                continue
                            if 'N' in don_mot_seq.upper():
                                continue

                            if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                                continue 

                            fseq_don = SeqRecord(don_mot_seq.upper(), id=fid, description='-ve label')
                            don_min_fh.write(fseq_don.format("fasta"))
                            break

                    elif loc[-1]=="-":
                        ## 
                        don_t_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]
                        idx = [xq.start() for xq in re.finditer(re.escape('AC'), str(don_t_seq).upper())]
                    
                        if boundary-1 in idx:
                            idx.remove(boundary-1)
                        if not idx:
                            continue

                        for xp in idx:
                            rloc_min = (int(loc[0])-boundary)+xp
                            don_mot_seq = rec.seq[(rloc_min-boundary)-1:(rloc_min+boundary)-1]
                            don_mot_seq = don_mot_seq.reverse_complement()
                        
                            if len(don_mot_seq) != 2*boundary:
                                continue 
                            if not don_mot_seq:
                                continue
                            if 'N' in don_mot_seq.upper():
                                continue

                            if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                                continue 

                            fseq_don = SeqRecord(don_mot_seq.upper(), id=fid, description='-ve label')
                            don_min_fh.write(fseq_don.format("fasta"))
                            #print don_mot_seq[0:10]
                            #print don_mot_seq[9:11]
                            break
                        ## 
                        acc_t_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                        idx = [xq.start() for xq in re.finditer(re.escape('CT'), str(acc_t_seq).upper())]

                        if boundary-1 in idx:
                            idx.remove(boundary-1)
                        if not idx:
                            continue

                        for xk in idx:
                            rloc_pos = (int(loc[1])-boundary)+xk
                            acc_mot_seq = rec.seq[(rloc_pos-boundary)+2:(rloc_pos+boundary)+2]
                            acc_mot_seq = acc_mot_seq.reverse_complement()
                    
                            if len(acc_mot_seq) != 2*boundary:
                                continue
                            if not acc_mot_seq:
                                continue
                            if 'N' in acc_mot_seq.upper():
                                continue

                            if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                                continue 

                            fseq_acc = SeqRecord(acc_mot_seq.upper(), id=fid, description='-ve label')
                            acc_min_fh.write(fseq_acc.format("fasta"))
                            break ## single label 
            
    foh.close()
    don_min_fh.close()
    acc_min_fh.close()

def min_seq_fetch(fnam, Label, signal, boundary):
    """
    fetch the signal sequence from the genome sequence
    """

    #TODO with reference to splice signals the TIS and TSS minus signal should follow the consensus region with true negative site in the sequence. 
    foh = helper._open_file(fnam)
    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam) ## result to respective dir 
    out_min_fh = open(out_path + "/" + signal + "_sig_minus_label.fa", 'w')

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                #try: ## random coordinate from mRNA regions particulary trying to exclude the real ones.
                    for ndr in range(4):
                        rloc = random.randint(loc[0]+boundary,loc[0]-boundary)
                        motif_seq = rec.seq[rloc-boundary:rloc+boundary]

                        if len(motif_seq) != 2*boundary:
                            continue
                        if not motif_seq:
                            continue
                        if 'N' in motif_seq:
                            continue

                        fseq = SeqRecord(motif_seq, id=fid+'_'+str(ndr), description='-ve label')
                        out_min_fh.write(fseq.format("fasta"))
                #except:
                #    rloc = random.randint(loc[0],loc[1])
                #    motif_seq = rec.seq[rloc-boundary:rloc+boundary]
                #    if not motif_seq:
                #        continue
                #    if 'N' in motif_seq:
                #        continue
                #    if len(motif_seq)!=400: 
                #        continue
                #    fseq = SeqRecord(motif_seq, id=fid, description='-ve label')
                #    out_min_fh.write(fseq.format("fasta"))
    out_min_fh.close()
    foh.close()

def true_ss_seq_fetch(fnam, Label, boundary):
    """
    True splice signals 
    """
    foh = helper._open_file(fnam)

    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    don_pos_fh = open(out_path + "/don_sig_plus_label.fa", 'w')
    acc_pos_fh = open(out_path + "/acc_sig_plus_label.fa", 'w')

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lfeat in Label[rec.id]:
                for fid, loc in Lfeat.items():
                    if loc[-1] == '+': 
                        acc_mot_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]
                        if len(acc_mot_seq) != 2*boundary:
                            continue
                        if not acc_mot_seq:
                            continue
                        if 'N' in acc_mot_seq.upper():
                            continue
                        if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                            continue 

                        fseq_acc = SeqRecord(acc_mot_seq.upper(), id=fid, description='+ve label')
                        acc_pos_fh.write(fseq_acc.format("fasta"))

                        don_mot_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                        if len(don_mot_seq) != 2*boundary:
                            continue 
                        if not don_mot_seq:
                            continue
                        if 'N' in don_mot_seq.upper():
                            continue
                        if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                            continue 

                        fseq_don = SeqRecord(don_mot_seq.upper(), id=fid, description='+ve label')
                        don_pos_fh.write(fseq_don.format("fasta"))

                    elif loc[-1] == '-':
                        don_mot_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]
                        don_mot_seq = don_mot_seq.reverse_complement()

                        if len(don_mot_seq) != 2*boundary:
                            continue 
                        if not don_mot_seq:
                            continue
                        if 'N' in don_mot_seq.upper():
                            continue
                        if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                            continue 

                        fseq_don = SeqRecord(don_mot_seq.upper(), id=fid, description='+ve label')
                        don_pos_fh.write(fseq_don.format("fasta"))

                        acc_mot_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                        acc_mot_seq = acc_mot_seq.reverse_complement()

                        if len(acc_mot_seq) != 2*boundary:
                            continue
                        if not acc_mot_seq:
                            continue
                        if 'N' in acc_mot_seq.upper():
                            continue
                        if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                            continue 

                        fseq_acc = SeqRecord(acc_mot_seq.upper(), id=fid, description='+ve label')
                        acc_pos_fh.write(fseq_acc.format("fasta"))

    don_pos_fh.close()
    acc_pos_fh.close()
    foh.close()

def pos_seq_fetch(fnam, Label, signal, boundary):
    """
    fetch the signal sequence from the genome sequence  
    """
    foh = helper._open_file(fnam)
    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    out_pos_fh = open(out_path + "/" + signal + "_sig_plus_label.fa", 'w')

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():

                    motif_seq = rec.seq[loc[0]-boundary:loc[0]+boundary]
                    if loc[-1] == '-': ## start codon position varies according to the strand
                        motif_seq = motif_seq.reverse_complement()

                    if len(motif_seq) != 2*boundary: 
                        continue
                    if not motif_seq:
                        continue
                    if 'N' in motif_seq:
                        continue

                    fseq = SeqRecord(motif_seq, id=fid, description='+ve label')
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
