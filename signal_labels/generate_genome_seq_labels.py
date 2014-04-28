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
import numpy 
import random
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from gfftools import helper, GFFParser 

def main(faname=None, gfname=None):
    """
    core unit

    @args faname: genome sequence in fasta 
    @type faname: str
    @args gfname: genome annotation in gff 
    @type gfname: str
    """
    
    # check for inputs
    if faname == gfname:
        print __doc__
        sys.exit(-1)

    # FIXME adjust the training label sequence count 
    label_cnt = 5000 # number of labels 

    # the number of positive and negative labels for training  
    plus_cnt = 1000
    minus_cnt = 10000

    # FIXME required input variables including the result path   
    #base_path = ''

    print 
    print 'processing genome annotation %s' % gfname
    print 
    anno_file_content = GFFParser.Parse(gfname) # extract genome annotation from gtf/gff type file 
    print 
    print 'processed genome annotation'
    print 
    
    # genomic signals : TranslationStop - TranscriptionStop - don/acc - Transcription - Translation 
    for signal in ['splice', 'cdsstop', 'cleave', 'tss', 'tis']: 

        gtf_db, feature_cnt, splice_checks, tid_gene_map = get_label_regions(anno_file_content, signal)
        print 
        print 'extracted %d %s signal regions' % (feature_cnt, signal)
        print 

        posLabel, COUNT = select_labels(gtf_db, feature_cnt, label_cnt) 
        print 
        print 'selecting %d RANDOM %s signal regions' % (COUNT, signal) 
        print 

        if signal == 'splice':
            acc_cnt, don_cnt = true_ss_seq_fetch(faname, posLabel) 
            print
            print 'selected %d acc %d don _plus_ %s signal lables' % (acc_cnt, don_cnt, signal)
            print

            acc_cnt, don_cnt = false_ss_seq_fetch(faname, posLabel, splice_checks, tid_gene_map)
            print 
            print 'selected %d acc %d don _minus_ %s signal lables' % (acc_cnt, don_cnt, signal)
            print 

        elif signal == 'tis':
            label_count = true_tis_seq_fetch(faname, posLabel)
            print 'selected %d plus %s signal lables' % (label_count, signal) 

            label_count = false_tis_seq_fetch(faname, posLabel)
            print 'selected %d minus %s signal lables' % (label_count, signal)

        elif signal == 'cdsstop':
            label_count = true_cdsStop_seq_fetch(faname, posLabel)
            print 'selected %d plus %s signal lables' % (label_count, signal)

            label_count = false_cdsStop_seq_fetch(faname, posLabel)
            print 'selected %d minus %s signal lables' % (label_count, signal)

        elif signal in ["tss", "cleave"]:
            label_count = plus_tss_seq_fetch(signal, faname, posLabel)
            print 'selected %d plus %s signal lables' % (label_count, signal) 

            label_count = minus_tss_seq_fetch(signal, faname, posLabel)
            print 'selected %d minus %s signal lables' % (label_count, signal) 

        # TODO remove the extra labels fetched from the previous step 
        plus_label_cleanup([signal], plus_cnt)

        minus_label_cleanup([signal], minus_cnt)

        # signal label processing over 
        print '%s signal done.' % signal
        print 
        break


def minus_label_cleanup(sig_type, minus_label_cnt):
    """
    clean up the result file with the correct number of fetched signal labels
    """

    #out_path = os.path.dirname(base_path) 
    sig_type = ['acc', 'don'] if sig_type[0] == "splice" else sig_type

    for signal in sig_type:

        # remove duplicate sequences, expecting the file to be in cwd path!  
        fh_seq = SeqIO.to_dict(SeqIO.parse("%s_sig_minus_label.fa" % signal, 'fasta')) 
        dup_ent = dict( (str(v.seq), k) for k,v in fh_seq.iteritems())
        non_dup_ent = dict((ele, 0) for ele in dup_ent.values())
        dup_ent.clear()

        #fasta_out_minus = open(out_path + "/"+ signal + "_sig_minus_label.bkp", 'w')
        fasta_out_minus = open("%s_sig_minus_label.bkp" % signal, 'w')
        cnt = 0 
    
        #minus_hd = open(out_path + "/" + signal +"_sig_minus_label.fa", "rU")
        minus_hd = open("%s_sig_minus_label.fa" % signal, "rU")
        for rec in SeqIO.parse(minus_hd, 'fasta'):
            if not rec.seq:
                continue
            # check for uniq sequence 
            if not rec.id in non_dup_ent:
                continue

            SeqIO.write([rec], fasta_out_minus, "fasta")
            cnt += 1
            if cnt == minus_label_cnt:
                break
        minus_hd.close()
        fasta_out_minus.close()
        # replacing with new file 
        #os.system('mv ' + out_path + '/'+ signal + '_sig_minus_label.bkp '+ out_path + "/" + signal + "_sig_minus_label.fa")
        os.system('mv %s_sig_minus_label.bkp %s_sig_minus_label.fa' % (signal, signal) )

        cnt = 0 
        # double checking the count  
        #minus_hd = open(out_path + "/" + signal + "_sig_minus_label.fa", "rU")
        minus_hd = open("%s_sig_minus_label.fa" % signal, "rU")
        for rec in SeqIO.parse(minus_hd, 'fasta'):
            cnt += 1
        minus_hd.close()

        print 'cleaned %d minus %s signal labels stored in %s_sig_minus_label.fa' % (cnt, signal, signal)


def plus_label_cleanup(sig_type, plus_label_cnt):
    """
    clean up the result file with the correct number of fetched signal labels
    """

    #out_path = os.path.dirname(base_path) 
    sig_type = ['acc', 'don'] if sig_type[0] == "splice" else sig_type

    for signal in sig_type:
        
        # remove duplicate sequences, expecting the file to be in cwd path!  
        fh_seq = SeqIO.to_dict(SeqIO.parse("%s_sig_plus_label.fa" % signal, 'fasta')) 
        dup_ent = dict( (str(v.seq), k) for k,v in fh_seq.iteritems())
        non_dup_ent = dict((ele, 0) for ele in dup_ent.values())
        dup_ent.clear()

        #fasta_out_plus = open(out_path + "/" + signal +"_sig_plus_label.bkp", 'w')
        fasta_out_plus = open("%s_sig_plus_label.bkp" % signal, 'w')

        cnt = 0 
        #plus_hd = open(out_path + "/"+ signal + "_sig_plus_label.fa", "rU")
        plus_hd = open("%s_sig_plus_label.fa" %signal, "rU")
        for rec in SeqIO.parse(plus_hd, 'fasta'):
            if not rec.seq:
                continue

            # check for uniq sequence 
            if not rec.id in non_dup_ent:
                continue

            SeqIO.write([rec], fasta_out_plus, "fasta")
            cnt += 1
            if cnt == plus_label_cnt:
                break
        plus_hd.close()
        fasta_out_plus.close()
        # replacing with new file 
        #os.system('mv ' + out_path + '/' + signal +'_sig_plus_label.bkp '+ out_path + "/"+ signal + "_sig_plus_label.fa")
        os.system('mv %s_sig_plus_label.bkp %s_sig_plus_label.fa' % (signal, signal) )
        
        cnt = 0 
        # double checking the count  
        #plus_hd = open(out_path + "/" + signal + "_sig_plus_label.fa", "rU")
        plus_hd = open("%s_sig_plus_label.fa" % signal, "rU")
        for rec in SeqIO.parse(plus_hd, 'fasta'):
            cnt += 1
        plus_hd.close()

        print 'cleaned %d plus %s signal labels stored in %s_sig_plus_label.fa' % (cnt, signal, signal)


def false_cdsStop_seq_fetch(fnam, Label, boundary=100):
    """
    fetch the minus cdsStop signal label sequences
    """
    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    #out_pos_fh = open(out_path + "/" + "tis_sig_plus_label.fa", 'w')
    out_min_fh = open("cdsstop_sig_minus_label.fa", 'w')

    true_label = 0 

    foh = helper._open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():

                    if loc[1] == '+': 
                        motif_seq = rec.seq[loc[2][0]:loc[2][1]]

                        # get index for negative signal label sequence site 
                        for ndr, stcodon in enumerate(['TAA', 'TAG', 'TGA']):
                            #search for different stop codon index 
                            idx = [xq.start() for xq in re.finditer(re.escape(stcodon), str(motif_seq).upper())]

                            if not idx:
                                continue

                            # get the false labels for randomly selected region or the available ones   
                            if len(idx) > 1:
                                idx = random.sample(idx, 1)

                            for nb, xp in enumerate(idx):

                                # adjusting the coordinate to the false site - mapping the relative position 
                                rloc_min = loc[2][0]+xp

                                # this has to make sure that we are removing the correct positive label 
                                if rloc_min == loc[0][0]:
                                    continue
                                
                                motif_sub_seq = rec.seq[(rloc_min-boundary)+1:(rloc_min+boundary)+2]

                                # check for sanity and consensus of the fetched sequence region 
                                if not motif_sub_seq:
                                    continue
                                if not all (XJ in 'ATCG' for XJ in str(motif_sub_seq.upper())):
                                    continue
                                if not str(motif_sub_seq[boundary-1:boundary+2]).upper() in ['TAA', 'TAG', 'TGA']:
                                    continue

                                # result to fasta out
                                fseq = SeqRecord(motif_sub_seq.upper(), id=fid+'_'+str(ndr)+'_'+str(nb), description='-ve label')
                                out_min_fh.write(fseq.format("fasta"))
                                true_label += 1 

                    elif loc[1] == '-': 
                        motif_seq = rec.seq[loc[2][0]:loc[2][1]]
                        
                        # 8713 3' UTR end, 3 nts for stop codon - 8714, 8715, 8716. 8717 end cds region
                        #print loc[0][0]-3 # 1-based coordinates 
                        #print loc[0][0] 
                        #print rec.seq[(loc[0][0]-4):(loc[0][0]-1)] # 0-based coordinates stop codon 

                        # get index for negative signal label sequence site 
                        for ndr, stcodon in enumerate(['TTA', 'TCA', 'CTA']):
                            #search for different stop codons 
                            idx = [xq.start() for xq in re.finditer(re.escape(stcodon), str(motif_seq).upper())]

                            if not idx:
                                continue
                                
                            # get the false labels for randomly selected region or the available ones   
                            if len(idx) > 1:
                                idx = random.sample(idx, 1)

                            for nb, xp in enumerate(idx):

                                # adjusting the coordinate to the false site 
                                rloc_min = loc[2][0]+xp
                                
                                #removing the true signal sequence site from selected false sites
                                if rloc_min == loc[0][0]-4:
                                    continue 

                                motif_sub_seq = rec.seq[(rloc_min-boundary)+1:(rloc_min+boundary)+2]
                                motif_sub_seq = motif_sub_seq.reverse_complement()

                                # check for sanity and consensus of the fetched sequence region 
                                if not motif_sub_seq:
                                    continue
                                if not all (XJ in 'ATCG' for XJ in str(motif_sub_seq.upper())):
                                    continue
                                if not str(motif_sub_seq[boundary-1:boundary+2]).upper() in ['TAA', 'TAG', 'TGA']:
                                    continue

                                # result to fasta out
                                fseq = SeqRecord(motif_sub_seq.upper(), id=fid+'_'+str(ndr)+'_'+str(nb), description='-ve label')
                                out_min_fh.write(fseq.format("fasta"))
                                true_label += 1 
    out_min_fh.close()
    foh.close()
    
    return true_label


def false_tis_seq_fetch(fnam, Label, boundary=100):
    """
    fetch the minus TIS signal label sequences 
    """

    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    #out_pos_fh = open(out_path + "/" + "tis_sig_plus_label.fa", 'w')
    out_min_fh = open("tis_sig_minus_label.fa", 'w')

    true_label = 0 

    foh = helper._open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():

                    if loc[1] == '+': 
                        motif_seq = rec.seq[loc[2][0]:loc[2][1]]

                        # get index for negative signal label sequence site 
                        idx = [xq.start() for xq in re.finditer(re.escape('ATG'), str(motif_seq).upper())]

                        if not idx:
                            continue 

                        # limit to take maximum 3 false labels from one defined feature
                        if len(idx) > 3:
                            idx = random.sample(idx, 3)

                        # get the false labels for randomly selected region or the available ones   
                        for ndr, xp in enumerate(idx):
                            # adjusting the coordinate to the false site 
                            rloc_min = loc[2][0]+xp

                            # removing the true signl sequence site from selected false sites
                            if rloc_min+1 == loc[0][0]:
                                continue

                            motif_seq = rec.seq[(rloc_min-boundary)+1:(rloc_min+boundary)+2]

                            # check for sanity and consensus of the fetched sequence region 
                            if len(motif_seq) != boundary*2+1:
                                continue
                            if not motif_seq:
                                continue
                            if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                                continue
                            if str(motif_seq[boundary-1:boundary+2]).upper() != 'ATG':
                                continue
                            # result to fasta out
                            fseq = SeqRecord(motif_seq.upper(), id=fid+'_'+str(ndr), description='-ve label')
                            out_min_fh.write(fseq.format("fasta"))
                            true_label += 1 

                    elif loc[1] == '-': 
                        motif_seq = rec.seq[loc[2][0]:loc[2][1]]

                        # get index for negative signal label sequence site 
                        idx = [xq.start() for xq in re.finditer(re.escape('CAT'), str(motif_seq).upper())]

                        if not idx:
                            continue 

                        # limit to take maximum 3 false labels from one defined feature
                        if len(idx) > 3:
                            idx = random.sample(idx, 3)

                        # get the false labels for randomly selected region or the available ones   
                        for ndr, xp in enumerate(idx):
                            # adjusting the coordinate to the false site 
                            rloc_min = loc[2][0]+xp

                            # removing the true signal sequence site from selected false sites
                            if rloc_min+3 == loc[0][0]:
                                continue
                            
                            motif_seq = rec.seq[(rloc_min-boundary)+1:(rloc_min+boundary)+2]
                            motif_seq = motif_seq.reverse_complement()

                            # check for sanity and consensus of the fetched sequence region 
                            if len(motif_seq) != boundary*2+1:
                                continue
                            if not motif_seq:
                                continue
                            if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                                continue
                            if str(motif_seq[boundary-1:boundary+2]).upper() != 'ATG':
                                continue

                            # result to fasta out
                            fseq = SeqRecord(motif_seq.upper(), id=fid+'_'+str(ndr), description='-ve label')
                            out_min_fh.write(fseq.format("fasta"))
                            true_label += 1 
    out_min_fh.close()
    foh.close()
    return true_label


def true_cdsStop_seq_fetch(fnam, Label, boundary=100):
    """
    fetch positive cdsStop signal labels 
    """

    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    #out_pos_fh = open(out_path + "/" + "cdsstop_sig_plus_label.fa", 'w')
    out_pos_fh = open("cdsstop_sig_plus_label.fa", 'w')

    true_label = 0 

    foh = helper._open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    if loc[1] == '+': 
                        motif_seq = rec.seq[(int(loc[0])-boundary)+1:int(loc[0])+boundary+2]

                        # check for sanity and consensus of the fetched sequence region 
                        if not motif_seq:
                            continue
                        if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                            continue
                        if not str(motif_seq[boundary-1:boundary+2]).upper() in ['TAA', 'TAG', 'TGA']:
                            continue

                        # result to fasta out
                        fseq = SeqRecord(motif_seq.upper(), id=fid, description='+ve label')
                        out_pos_fh.write(fseq.format("fasta"))
                        true_label += 1

                    elif loc[1] == '-': 
                        motif_seq = rec.seq[(int(loc[0])-boundary)-3:(int(loc[0])+boundary)-2]
                        motif_seq = motif_seq.reverse_complement()

                        # check for sanity and consensus of the fetched sequence region 
                        if not motif_seq:
                            continue
                        if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                            continue
                        if not str(motif_seq[boundary-1:boundary+2]).upper() in ['TAA', 'TAG', 'TGA']:
                            continue

                        # result to fasta out
                        fseq = SeqRecord(motif_seq.upper(), id=fid, description='+ve label')
                        out_pos_fh.write(fseq.format("fasta"))
                        true_label += 1
    out_pos_fh.close()
    foh.close()
    return true_label


def true_tis_seq_fetch(fnam, Label, boundary=100):
    """
    fetch the plus TIS signal sequence.
    """

    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    #out_pos_fh = open(out_path + "/" + "tis_sig_plus_label.fa", 'w')
    out_pos_fh = open("tis_sig_plus_label.fa", 'w')

    true_label = 0 

    foh = helper._open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():

                    if loc[1] == '+': 
                        motif_seq = rec.seq[int(loc[0])-boundary:int(loc[0])+boundary+1]

                        # check for sanity and consensus of the fetched sequence region 
                        if len(motif_seq) != boundary*2+1:
                            continue
                        if not motif_seq:
                            continue
                        if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                            continue
                        if str(motif_seq[boundary-1:boundary+2]).upper() != 'ATG':
                            continue

                        # result to fasta out
                        fseq = SeqRecord(motif_seq.upper(), id=fid, description='+ve label')
                        out_pos_fh.write(fseq.format("fasta"))
                        true_label += 1 

                    elif loc[1] == '-': 
                        motif_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-1]
                        motif_seq = motif_seq.reverse_complement()

                        # check for sanity and consensus of the fetched sequence region 
                        if len(motif_seq) != boundary*2+1:
                            continue
                        if not motif_seq:
                            continue
                        if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                            continue
                        if str(motif_seq[boundary-1:boundary+2]).upper() != 'ATG':
                            continue

                        # result to fasta out
                        fseq = SeqRecord(motif_seq.upper(), id=fid, description='+ve label')
                        out_pos_fh.write(fseq.format("fasta"))
                        true_label += 1 
    out_pos_fh.close()
    foh.close()
    return true_label 

     
def get_label_regions(gtf_content, signal):
    """
    get signal sequence location from the annotation
    
    @args gtf_content: parsed object from gtf/gff file 
    @type gtf_content: numpy array 
    @args signal: genomic signal name 
    @type signal: str 
    """

    feat_cnt = 0
    anno_db = defaultdict(list) 
    splice_signal_point = defaultdict(list)
    trans_gene_map = dict() 
    
    for feature in gtf_content: # gene list from GFFParse function 
        mod_anno_db = dict()

        if signal == 'tss':
            for xp, ftid in enumerate(feature['transcripts']):
                feat_cnt += 1
                mod_anno_db[ftid[0]] = (feature['tss'][xp], 
                                feature['strand'], 
                                (int(feature['start']), int(feature['stop']))
                                )
        elif signal == 'tis':
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
        elif signal == 'cleave':
            for xp, ftid in enumerate(feature['transcripts']):
                feat_cnt += 1
                mod_anno_db[ftid[0]] = (feature['cleave'][xp], 
                                feature['strand'],
                                (int(feature['start']), int(feature['stop']))
                                )
        elif signal == 'splice':
            for xp, ftid in enumerate(feature['transcripts']):# considering each annotated transcript 
                exons_count = len(feature['exons'][xp])
                
                if exons_count > 2:# spliced transcripts 
                    id_cnt = 1
                    exon_ends = [] 
                    for ikx, ex in enumerate(feature['exons'][xp]):
                        feat_cnt += 1
                        exon_ends.extend((ex[0], ex[1]))

                        alt_id = "%s.%d" % (ftid[0], id_cnt)
                        id_cnt += 1
                        
                        if ikx == 0: # first and last exon ends are not considering 
                            ex[0] = None
                        if exons_count-1 == ikx: 
                            ex[1] = None 
                        
                        mod_anno_db[alt_id] = (ex[0], 
                                            ex[1], 
                                            feature['strand']
                                            )
                    splice_signal_point[feature['name']].extend(exon_ends) # appending splice sites of transcripts to the gene
                    trans_gene_map[ftid[0]] = feature['name'] # transcript to gene mapping 

        if mod_anno_db:
            anno_db[feature['chr']].append(mod_anno_db)
    
    return dict(anno_db), feat_cnt, splice_signal_point, trans_gene_map 


def false_ss_seq_fetch(fnam, Label, don_acc_check, tr_gene_mp, boundary=100, sample=10):
    """
    false splice signals
    """

    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam) 
    #don_min_fh = open(out_path + "/don_sig_minus_label.fa", 'w')
    #acc_min_fh = open(out_path + "/acc_sig_minus_label.fa", 'w')
    don_min_fh = open("don_sig_minus_label.fa", 'w')
    acc_min_fh = open("acc_sig_minus_label.fa", 'w')

    true_label_acc = 0 
    true_label_don = 0 

    foh = helper._open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for fLabel in Label[rec.id]:
                for fid, loc in fLabel.items():
                    
                    fid_lookup = re.search(r'^(.*)\.\d+$', fid).group(1) # to search in tr_gene_map 
                    signal_location = don_acc_check[tr_gene_mp[fid_lookup]]
                    signal_location.sort() 

                    if loc[-1]=='+':
                        if not numpy.isnan(loc[0]):
                            # acceptor splice site signal 
                            acc_t_seq = rec.seq[int(signal_location[0]):int(signal_location[-1])] # from gene body 
                            # get index for negative signal label sequence site 
                            idx = [xq.start() for xq in re.finditer(re.escape('AG'), str(acc_t_seq).upper())]

                            # limit to take maximum 2 false labels from one defined feature
                            if len(idx) > sample:
                                idx = random.sample(idx, sample)
                            
                            for ndr, xp in enumerate(idx):
                                # adjusting the coordinate to the false site 
                                rloc_min = int(signal_location[0])+xp
                                if rloc_min+3 in signal_location: # ---AG|Exon~~~~ removing the true signal from false site 
                                    continue

                                acc_mot_seq = rec.seq[(rloc_min-boundary)+1:(rloc_min+boundary)+1]

                                # check for sanity and consensus of the fetched sequence region 
                                if len(acc_mot_seq) != 2*boundary:
                                    continue
                                if not acc_mot_seq:
                                    continue
                                if not all (XJ in 'ATCG' for XJ in str(acc_mot_seq.upper())):
                                    continue
                                if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                                    continue 

                                # write to fasta out 
                                fseq_acc = SeqRecord(acc_mot_seq.upper(), id='%s_%d' % (fid, ndr), description='-ve label')
                                acc_min_fh.write(fseq_acc.format("fasta"))
                                true_label_acc += 1 
                        
                        if not numpy.isnan(loc[1]):
                            # donor splice site signal 
                            don_t_seq = rec.seq[int(signal_location[0]):int(signal_location[-1])] # from gene body 
                            # get index for negative signal label sequence site 
                            idx = [xq.start() for xq in re.finditer(re.escape('GT'), str(don_t_seq).upper())]

                            # limit to take maximum 2 false labels from one defined feature
                            if len(idx) > sample:
                                idx = random.sample(idx, sample)


                            for ndr, xj in enumerate(idx):
                                # adjusting the coordinate to the false site 
                                rloc_pos = int(signal_location[0])+xj
                                if rloc_pos in signal_location:
                                    continue

                                don_mot_seq = rec.seq[(rloc_pos-boundary)+1:(rloc_pos+boundary)+1]

                                # check for sanity and consensus of the fetched sequence region 
                                if len(don_mot_seq) != 2*boundary:
                                    continue
                                if not don_mot_seq:
                                    continue
                                if not all (XJ in 'ATCG' for XJ in str(don_mot_seq.upper())):
                                    continue
                                if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                                    continue 

                                # write to fasta out 
                                fseq_don = SeqRecord(don_mot_seq.upper(), id='%s_%d' % (fid, ndr), description='-ve label')
                                don_min_fh.write(fseq_don.format("fasta"))
                                true_label_don += 1 

                    elif loc[-1]=="-":
                        if not numpy.isnan(loc[0]):
                            # donor splice signal site 
                            don_t_seq = rec.seq[int(signal_location[0]):int(signal_location[-1])]
                            # get index for negative signal label sequence site 
                            idx = [xq.start() for xq in re.finditer(re.escape('AC'), str(don_t_seq).upper())]

                            # limit to take maximum 2 false labels from one defined feature
                            if len(idx) > sample:
                                idx = random.sample(idx, sample)

                            for ndr, xp in enumerate(idx):
                                # adjusting the coordinate to the false site 
                                rloc_min = int(signal_location[0])+xp
                                if rloc_min+3 in signal_location:
                                    continue

                                don_mot_seq = rec.seq[(rloc_min-boundary)+1:(rloc_min+boundary)+1]
                                don_mot_seq = don_mot_seq.reverse_complement()

                                # check for sanity and consensus of the fetched sequence region 
                                if len(don_mot_seq) != 2*boundary:
                                    continue 
                                if not don_mot_seq:
                                    continue
                                if not all (XJ in 'ATCG' for XJ in str(don_mot_seq.upper())):
                                    continue
                                if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                                    continue 

                                # write to fasta out 
                                fseq_don = SeqRecord(don_mot_seq.upper(), id='%s_%d' % (fid, ndr), description='-ve label')
                                don_min_fh.write(fseq_don.format("fasta"))
                                true_label_don += 1 
                        
                        if not numpy.isnan(loc[1]):
                            # acceptor splice signal site 
                            acc_t_seq = rec.seq[int(signal_location[0]):int(signal_location[-1])]
                            # get index for negative signal label sequence site 
                            idx = [xq.start() for xq in re.finditer(re.escape('CT'), str(acc_t_seq).upper())]

                            # limit to take maximum 2 false labels from one defined feature
                            if len(idx) > sample:
                                idx = random.sample(idx, sample)

                            for ndr, xk in enumerate(idx):
                                # adjusting the coordinate to the false site 
                                rloc_pos = int(signal_location[0])+xk
                                if rloc_pos in signal_location:
                                    continue

                                acc_mot_seq = rec.seq[(rloc_pos-boundary)+1:(rloc_pos+boundary)+1]
                                acc_mot_seq = acc_mot_seq.reverse_complement()

                                # check for sanity and consensus of the fetched sequence region 
                                if len(acc_mot_seq) != 2*boundary:
                                    continue
                                if not acc_mot_seq:
                                    continue
                                if not all (XJ in 'ATCG' for XJ in str(acc_mot_seq.upper())):
                                    continue
                                if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                                    continue 

                                # write to fasta out 
                                fseq_acc = SeqRecord(acc_mot_seq.upper(), id='%s_%d' % (fid, ndr), description='-ve label')
                                acc_min_fh.write(fseq_acc.format("fasta"))
                                true_label_acc += 1

    foh.close()
    don_min_fh.close()
    acc_min_fh.close()
    return true_label_acc, true_label_don


def minus_tss_seq_fetch(signal, fnam, Label, boundary=100):
    """
    fetch the minus TSS, cleave signal sequence label
    """

    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam) ## result to fasta base dir 
    #out_min_fh = open(out_path + "/" + "_sig_minus_label.fa", 'w')
    out_min_fh = open(signal + "_sig_minus_label.fa", 'w')

    true_label = 0 

    foh = helper._open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    for ndr in range(3):

                        # max 3 labels from a true feature 
                        rloc = random.randint(loc[2][0],loc[2][1])

                        # remove the true signal index from random sampling 
                        if rloc == int(loc[0]):
                            continue

                        motif_seq = rec.seq[rloc-boundary:rloc+boundary+1]

                        # sanity check for featched sequence 
                        if len(motif_seq) != 2*boundary+1:
                            continue
                        if not motif_seq:
                            continue
                        if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                            continue

                        # write to fasta out 
                        fseq = SeqRecord(motif_seq.upper(), id=fid+'_'+str(ndr), description='-ve label')
                        out_min_fh.write(fseq.format("fasta"))
                        true_label += 1 
    out_min_fh.close()
    foh.close()
    return true_label


def true_ss_seq_fetch(fnam, Label, boundary=100):
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

    true_label_acc = 0
    true_label_don = 0 

    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lfeat in Label[rec.id]:
                for fid, loc in Lfeat.items():
                    if loc[-1] == '+': 
                        if not numpy.isnan(loc[0]):
                            # acceptor splice site 
                            acc_mot_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]

                            # sanity check and consensus sequence from fetched sequence 
                            if len(acc_mot_seq) != boundary*2:
                                continue
                            if not acc_mot_seq:
                                continue
                            if not all (XJ in 'ATCG' for XJ in str(acc_mot_seq.upper())):
                                continue
                            if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                                continue 

                            # write to fasta out 
                            fseq_acc = SeqRecord(acc_mot_seq.upper(), id=fid, description='+ve label')
                            acc_pos_fh.write(fseq_acc.format("fasta"))
                            true_label_acc += 1 
                        
                        if not numpy.isnan(loc[1]):
                            # donor splice site 
                            don_mot_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]

                            if len(don_mot_seq) != 2*boundary:
                                continue 
                            if not don_mot_seq:
                                continue
                            if not all (XJ in 'ATCG' for XJ in str(don_mot_seq.upper())):
                                continue
                            if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                                continue 

                            # write to fasta out 
                            fseq_don = SeqRecord(don_mot_seq.upper(), id=fid, description='+ve label')
                            don_pos_fh.write(fseq_don.format("fasta"))
                            true_label_don += 1

                    elif loc[-1] == '-':
                        if not numpy.isnan(loc[0]):
                            # donor splice site signal 
                            don_mot_seq = rec.seq[(int(loc[0])-boundary)-2:(int(loc[0])+boundary)-2]
                            don_mot_seq = don_mot_seq.reverse_complement()

                            # sanity check and consensus sequence from fetched sequence 
                            if len(don_mot_seq) != 2*boundary:
                                continue 
                            if not don_mot_seq:
                                continue
                            if not all (XJ in 'ATCG' for XJ in str(don_mot_seq.upper())):
                                continue
                            if str(don_mot_seq[boundary-1:boundary+1]).upper() != 'GT':
                                continue 

                            # write to fasta out 
                            fseq_don = SeqRecord(don_mot_seq.upper(), id=fid, description='+ve label')
                            don_pos_fh.write(fseq_don.format("fasta"))
                            true_label_don += 1 
                        
                        if not numpy.isnan(loc[1]):
                            # acceptor splice signal 
                            acc_mot_seq = rec.seq[(int(loc[1])-boundary)+1:(int(loc[1])+boundary)+1]
                            acc_mot_seq = acc_mot_seq.reverse_complement()

                            # sanity check and consensus sequence from fetched sequence 
                            if len(acc_mot_seq) != 2*boundary:
                                continue
                            if not acc_mot_seq:
                                continue
                            if not all (XJ in 'ATCG' for XJ in str(acc_mot_seq.upper())):
                                continue
                            if str(acc_mot_seq[boundary-1:boundary+1]).upper() != 'AG':
                                continue 

                            # write to fasta out 
                            fseq_acc = SeqRecord(acc_mot_seq.upper(), id=fid, description='+ve label')
                            acc_pos_fh.write(fseq_acc.format("fasta"))
                            true_label_acc += 1 
    don_pos_fh.close()
    acc_pos_fh.close()
    foh.close()
    return true_label_acc, true_label_don


def plus_tss_seq_fetch(signal, fnam, Label, boundary=100):
    """
    fetch the plus TSS, cleave signal sequence. The default flanking region is 200 nucleotides. 
    """

    real_fnam = os.path.realpath(fnam)
    out_path = os.path.dirname(real_fnam)
    #out_pos_fh = open(out_path + "/" + "_sig_plus_label.fa", 'w')
    out_pos_fh = open(signal + "_sig_plus_label.fa", 'w')

    true_label = 0 

    foh = helper._open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    motif_seq = rec.seq[int(loc[0])-boundary:int(loc[0])+boundary+1]

                    if loc[1] == '-': 
                        motif_seq = motif_seq.reverse_complement()

                    # sanity check for the fetched sequence 
                    if len(motif_seq) != boundary*2+1: 
                        continue
                    if not motif_seq:
                        continue
                    if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                        continue

                    # write to fasta out 
                    fseq = SeqRecord(motif_seq.upper(), id=fid, description='+ve label')
                    out_pos_fh.write(fseq.format("fasta"))
                    true_label += 1 

    out_pos_fh.close()
    foh.close()
    return true_label


def select_labels(feat_db, feat_count, label_cnt):
    """
    Random sampling to select signal lables
    """

    assert label_cnt <= feat_count, 'Number of features annotated %d' % feat_count

    try:
        accept_prob = (1.0*label_cnt)/feat_count
    except:
        accept_prob = 1

    while 1: # ensure the label count 
        counter, LSet = recursive_fn(feat_db, label_cnt, accept_prob)
        if label_cnt <= counter:
            break
        print '    recursive ... %d' % counter

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

    try:
        fas_name = sys.argv[1]
        gff_name = sys.argv[2]
    except:
        print __doc__
        sys.exit(-1) 

    main(fas_name, gff_name)
