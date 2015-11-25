#!/usr/bin/env python
"""
Extract genomic signal label sequences from the genome sequence and annotation file. 
    Transcription Start Site[tss], 
    Translational Initiation Site[tis], 
    Splice Sites [don, acc], 
    Transcription Stop Site[cleave]
    Translation Stop Site [cdsstop]

    the resulting label sequences are placed under the working dir in fasta format.

The signal sequence are flanked by defined number of nucleotides upstream and downstream of genome signal region.

Usage: 
    python fetch_signal_seq_examples.py in.fasta.(gz|bz2) in.gtf.(gz|bz2)

Requirements:
    BioPython:- http://biopython.org 
    gfftools:- https://github.com/vipints/gfftools
"""

import re
import os 
import sys
import numpy 
import shutil
import random
import filecmp
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from gfftools import helper, GFFParser 


def get_feature_seq(faname, gfname, signal='tss', label_cnt=4000, plus_cnt=3000, minus_cnt=9000, flanks=1200):
    """
    core unit to fetch genomic signal label sequences 

    @args faname: genome sequence in fasta 
    @type faname: str
    @args gfname: genome annotation in gff 
    @type gfname: str
    @args signal: different genomic signal types ['splice', 'cdsstop', 'cleave', 'tss', 'tis']
    @type signal: str
    @args label_cnt: number of random labels to be selected 
    @type label_cnt: integer
    @args plus_cnt: plus label count TODO default parameter details 
    @type plus_cnt: integer
    @args minus_cnt: minus label count 
    @type minus_cnt: integer
    @args flanks: flanking sequence length 
    @type flanks: integer
    """
    
    ## check the files
    if filecmp.cmp(faname, gfname):
        exit("Do the two files are exactly same? Please check that!")

    ##FIXME required input variables including the result path   
    #base_path = ''

    ## extract genome annotation from gtf/gff type file 
    sys.stdout.write('reading genome annotation from %s...' % gfname)
    anno_file_content = GFFParser.Parse(gfname) 
    sys.stdout.write(' ...done.\n')

    ## check the consistency of chr names in fasta and gff file 
    chrom_name_consistency(faname, anno_file_content) 
    
    gtf_db, feature_cnt, signal_checks, tid_gene_map = get_label_regions(anno_file_content, signal)
    sys.stdout.write('extracted %d %s signal regions\n' % (feature_cnt, signal))

    posLabel, COUNT = select_labels(gtf_db, feature_cnt, label_cnt) 
    sys.stdout.write('selecting %d RANDOM %s signal regions\n' % (COUNT, signal))

    if signal == "tss": 
        label_count_plus = plus_tss_cleave_seq_fetch(signal, faname, posLabel, flanks)
        print 'selected %d positive %s signal examples' % (label_count_plus, signal) 

        label_count_minus = minus_tss_seq_fetch(faname, posLabel, signal_checks, tid_gene_map, flanks)
        print 'selected %d negative %s signal examples' % (label_count_minus, signal) 

        #FIXME cleanup 
        label_ids_plus = plus_label_cleanup([signal], plus_cnt, label_count_plus)
        label_id_minus = minus_label_cleanup([signal], minus_cnt, label_count_minus)

        ## check for the consistent examples ratio from plus and minus set  
        label_count_pos = strip_examples_without_proper_pairs(signal) 

        ## remove the extra labels fetched from the previous step 
        matching_pos_neg_example(signal, label_count_pos, plus_cnt)
         
    # signal label processing over 
    print '%s signal done.' % signal
    print 


def strip_examples_without_proper_pairs(signal):
    """
    This helps to get the final list of examples with 1:1 ratio of pos and neg examples 
    """

    min_pos_neg_ratio = 3 

    in_fas_pos = "%s_sig_pos_example.fa" % signal
    in_fas_neg = "%s_sig_neg_example.fa" % signal

    pos_fas_rec = defaultdict(list) 
    for rec in SeqIO.parse(in_fas_pos, 'fasta'):
        desc = rec.description.split(" ")
        pos_fas_rec[desc[-1]].append(rec.description) ## FIXME if the transcript_id is not uniq 

    pos_neg_ratio = defaultdict(int) 
    for rec_neg in SeqIO.parse(in_fas_neg, 'fasta'):
        desc_neg = rec_neg.description.split(" ")
        if desc_neg[-1] in pos_fas_rec:
            pos_neg_ratio[desc_neg[-1]] += 1
    
    ## change the pos examples and write a new file 
    pos_rec_fas_file = "%s_sig_pos_example.bkp" % signal
    pos_rec_fas_fh = open(pos_rec_fas_file, "w") 
    pair_pos_example = 0 
    for rec in SeqIO.parse(in_fas_pos, 'fasta'):
        desc = rec.description.split(" ")
        if desc[-1] in pos_neg_ratio:
            if pos_neg_ratio[desc[-1]] != min_pos_neg_ratio:
                continue
            pos_rec_fas_fh.write(rec.format("fasta"))
            pair_pos_example += 1 
             
    pos_rec_fas_fh.close()
    print "%d positive examples with proper pairs" % pair_pos_example

    shutil.move('%s_sig_pos_example.bkp' % signal, '%s_sig_pos_example.fa' % signal)

    return pair_pos_example


def iterative_fetch(in_fas_pos, out_fas_pos, sub_sample_records, accept_prob):
    """
    the iterative mode of selecting records to match the required number

    @args in_fas_pos:
    @type in_fas_pos:
    @args out_fas_pos:
    @type out_fas_pos:
    @args sub_sample_records:
    @type sub_sample_records:
    @args accept_prob:
    @type accept_prob:
    """

    fasta_rec = 0 
    sample_cnt = 1
    pos_fas_rec = defaultdict(list) 
    pos_fas_out = open(out_fas_pos, "w") 

    for rec in SeqIO.parse(in_fas_pos, 'fasta'):
        fasta_rec += 1 
        desc = rec.description.split(" ")

        rnb = random.random()
        if rnb <= accept_prob:
            pos_fas_out.write(rec.format("fasta"))
            pos_fas_rec[desc[-1]].append(rec.description) ## transcript_id is uniq. 
            if sample_cnt == sub_sample_records:
                break 
            sample_cnt += 1 

    pos_fas_out.close()

    return fasta_rec, sample_cnt, pos_fas_rec
   

def matching_pos_neg_example(signal, total_record_count, sub_sample_records):
    """
    fetch random pos examples and corresponding neg example in 1:1 ratio 
    ## based on the pos examples 
    """

    ## FIXME this has to be automatized 
    neg_sample_ratio = 3 ## sub sample ratio pos:neg  

    in_fas_pos = "%s_sig_pos_example.fa" % signal
    in_fas_neg = "%s_sig_neg_example.fa" % signal

    try:
        accept_prob = (1.0*sub_sample_records)/total_record_count
    except:
        accept_prob = 1

    out_fas_pos = "%s_sig_pos_example.bkp" % signal

    while True:
        fasta_rec, sample_cnt, pos_fas_rec = iterative_fetch(in_fas_pos, out_fas_pos, sub_sample_records, accept_prob) 
        if sample_cnt >= sub_sample_records:
            break
        sys.stdout.write('still trying ... %d\n' % sample_cnt) 
 
    sys.stdout.write('%d number of records scanned\n' % fasta_rec)
    sys.stdout.write('%d number of pos examples dumped\n' % sample_cnt)

    shutil.move('%s_sig_pos_example.bkp' % signal, '%s_sig_pos_example.fa' % signal)

    neg_rec_cnt = 0 
    pos_neg_ratio = defaultdict(int) 
    out_fas_neg = "%s_sig_neg_example.bkp" % signal 
    neg_fas_out = open(out_fas_neg, "w") 
    for rec_neg in SeqIO.parse(in_fas_neg, 'fasta'):
        desc_neg = rec_neg.description.split(" ")

        if desc_neg[-1] in pos_fas_rec:
            pos_neg_ratio[desc_neg[-1]] += 1
             
            if pos_neg_ratio[desc_neg[-1]] > neg_sample_ratio:
                continue 

            neg_fas_out.write(rec_neg.format("fasta"))
            neg_rec_cnt += 1 

    neg_fas_out.close() 
    sys.stdout.write('%d number of neg examples dumped\n' % neg_rec_cnt)

    shutil.move('%s_sig_neg_example.bkp' % signal, '%s_sig_neg_example.fa' % signal)


def chrom_name_consistency(fasta_fname, gff_content):
    """
    check the chromosome name consistency in gff file and fasta file 

    @args fasta_fname: fasta file 
    @type fasta_fname: string 
    @args gff_content: parsed object from gtf/gff file 
    @type gff_content: numpy array 
    """
    
    chrom_fasta = dict() 
    fasta_h = helper.open_file(fasta_fname)
    for fas_rec in SeqIO.parse(fasta_h, "fasta"):
        chrom_fasta[fas_rec.id] = 0 

    chrom_gff = dict() 
    for gff_rec in gff_content:
        chrom_gff[gff_rec['chr']] = 0 

    cnt = 0 
    for chrom in chrom_fasta:
        if not chrom in chrom_gff:
            sys.stdout.write('%s NOT found in GFF/GTF file\n' % chrom)
            cnt += 1  
    
    if cnt == len(chrom_fasta):
        exit('error: chromosome/contig names are different in provided fasta, gff/gtf file.')


def fetch_unique_labels(firstLabel, secondLabel):
    """
    fetch separate set of negative training labels 

    @args firstLabel: identifiers from positive labels. example {ATG001:0, ATG002:0}
    @type firstLabel: dict() 
    @args secondLabel: random labels selected for the second time. example chr2:[{ATG001:(10, 20, +), ATG002:(20, 40, -)}]
    @type secondLabel: defaultdict(list) 
    """
    unique_label_loc = defaultdict(list) 

    for chrom, feat_loc_2 in secondLabel.items():
        #TODO better way, getting the uniq ids from second round of labels 
        diff_keys = set(feat_loc_2[0].keys()) - set(firstLabel.keys()) 
        ## fetch those keys details
        ftdet = dict() 
        for ft in list(diff_keys):
            if ft in feat_loc_2[0]:
                ftdet[ft] = feat_loc_2[0][ft] 
        unique_label_loc[chrom].append(ftdet) 
           
    return unique_label_loc



def minus_label_cleanup(sig_type, minus_label_cnt, feat_count):
    """
    clean up the result file with the correct number of fetched signal labels

    @args sig_type: genomic signal type 
    @type sig_type: str 
    @args plus_label_cnt: number of positive labels 
    @type plus_label_cnt: int 
    @args feat_count: total number of random valid features extracted
    @type feat_count: int 
    """

    assert minus_label_cnt <= feat_count, 'PLEASE INCREASE THE NUMBER OF RANDOM FEATURES TO BE SELECTED FROM %d' % feat_count

    #out_path = os.path.dirname(base_path) 
    sig_type = ['acc', 'don'] if sig_type[0] == "splice" else sig_type

    for signal in sig_type:
        # 1 remove the same rec.id presents in the file 
        fas_rec_ids = dict() 
        for fas_rec in SeqIO.parse("%s_sig_neg_example.fa" % signal, 'fasta'):
            fas_rec_ids[fas_rec.id] = 0 
        
        # remove the duplicated records based on the record name 
        fasta_out = open("%s_sig_neg_example.bkp" % signal, 'w')
        for fas_rec in SeqIO.parse("%s_sig_neg_example.fa" % signal, 'fasta'):
            if fas_rec.id in fas_rec_ids:
                SeqIO.write([fas_rec], fasta_out, "fasta")
                del fas_rec_ids[fas_rec.id]
        fasta_out.close() 
        shutil.move('%s_sig_neg_example.bkp' % signal, '%s_sig_neg_example.fa' % signal)

        # 2 remove duplicate sequences, expecting the file to be in cwd path!  
        fh_seq = SeqIO.to_dict(SeqIO.parse("%s_sig_neg_example.fa" % signal, 'fasta')) 
        dup_ent = dict( (str(v.seq), k) for k,v in fh_seq.iteritems())
        non_dup_ent = dict((ele, 0) for ele in dup_ent.values())
        dup_ent.clear()
        
        assert minus_label_cnt < len(non_dup_ent), 'DUPLICATE ENTRIES PRESENT NON-DUPLICATE ONES ARE %d' % len(non_dup_ent)  
        counter, label_seq_ids = random_pick(signal, 'neg', non_dup_ent, minus_label_cnt)
        
        #os.system('mv ' + out_path + '/'+ signal + '_sig_minus_label.bkp '+ out_path + "/" + signal + "_sig_minus_label.fa")
        shutil.move('%s_sig_neg_example.bkp' % signal, '%s_sig_neg_example.fa' % signal )
        sys.stdout.write('cleaned %d minus %s signal examples stored in %s_sig_neg_example.fa\n' % (counter, signal, signal))

        return label_seq_ids


def plus_label_cleanup(sig_type, plus_label_cnt, feat_count):
    """
    clean up the result file with the correct number of fetched signal labels
 
    @args sig_type: genomic signal type 
    @type sig_type: str 
    @args plus_label_cnt: number of positive labels 
    @type plus_label_cnt: int 
    @args feat_count: total number of random valid features extracted
    @type feat_count: int 
    """
    assert plus_label_cnt <= feat_count, 'PLEASE INCREASE THE NUMBER OF RANDOM FEATURES TO BE SELECTED FROM %d' % feat_count

    sig_type = ['acc', 'don'] if sig_type[0] == "splice" else sig_type
    for signal in sig_type:
        # 1 remove the same rec.id presents in the file 
        fas_rec_ids = dict() 
        for fas_rec in SeqIO.parse("%s_sig_pos_example.fa" % signal, 'fasta'):
            fas_rec_ids[fas_rec.id] = 0 
        
        # remove the duplicated records based on the record name 
        fasta_out = open("%s_sig_pos_example.bkp" % signal, 'w')
        for fas_rec in SeqIO.parse("%s_sig_pos_example.fa" % signal, 'fasta'):
            if fas_rec.id in fas_rec_ids:
                SeqIO.write([fas_rec], fasta_out, "fasta")
                del fas_rec_ids[fas_rec.id]
        fasta_out.close() 
        shutil.move('%s_sig_pos_example.bkp' % signal, '%s_sig_pos_example.fa' % signal)

        # 2 remove duplicate sequences, expecting the file to be in cwd path!  
        fh_seq = SeqIO.to_dict(SeqIO.parse("%s_sig_pos_example.fa" % signal, 'fasta')) 
        dup_ent = dict((str(v.seq), k) for k,v in fh_seq.iteritems())
        non_dup_ent = dict((ele, 0) for ele in dup_ent.values())
        dup_ent.clear()

        assert plus_label_cnt < len(non_dup_ent), 'DUPLICATE ENTRIES PRESENT NON-DUPLICATE ONES ARE %d' % len(non_dup_ent)  
        #FIXME cleanup 
        counter, label_seq_ids = random_pick(signal, 'pos', non_dup_ent, plus_label_cnt)

        shutil.move('%s_sig_pos_example.bkp' % signal, '%s_sig_pos_example.fa' % signal)
        sys.stdout.write('cleaned %d plus %s signal examples stored in %s_sig_pos_example.fa\n' % (counter, signal, signal))

        return label_seq_ids


def random_pick(signal, plus_minus, non_dup_ent, lb_cnt):
    """
    navigate through the filtered set and fetch the right number of plus and minus labels 

    @args signal: genomic signal type  
    @type signal: str 
    @args plus_minus: plus or minus tag 
    @type plus_minus: str 
    @args non_dup_ent: unique entries 
    @type non_dup_ent: dict 
    @args lb_cnt: plus or minus label count  
    @type lb_cnt: int  
    """

    cnt = 0 
    label_seq_ids = dict()

    #fasta_out_plus = open(out_path + "/" + signal +"_sig_plus_label.bkp", 'w')
    fasta_out_plus = open("%s_sig_%s_example.bkp" % (signal, plus_minus), 'w')

    #plus_hd = open(out_path + "/"+ signal + "_sig_plus_label.fa", "rU")
    plus_hd = open("%s_sig_%s_example.fa" % (signal, plus_minus), "rU")
    for rec in SeqIO.parse(plus_hd, 'fasta'):
        if not rec.seq:
            continue

        # check for uniq sequence 
        if not rec.id in non_dup_ent:
            continue

        SeqIO.write([rec], fasta_out_plus, "fasta")
        cnt += 1

        desc = rec.description.split(' ')
        #FIXME cleanup 
        label_seq_ids[desc[-1]] = 0

    plus_hd.close()
    fasta_out_plus.close()
    
    return cnt, label_seq_ids 


def false_cdsStop_seq_fetch(fnam, Label, cdsstop_check, tr_gene_mp, boundary=100, sample=2):
    """
    fetch the minus cdsStop signal label sequences

    @args fnam: genome sequence in fasta format 
    @type fnam: str
    @args Label: signal sequence in the genome <chrom<transcript:(location, strand)>>
    @type Label: dict 
    @args cdsstop_check: cdsstop signals from different transcripts of the same gene <transcript:[exon_end_locations]>  
    @type cdsstop_check: defaultdict(list)
    @args tr_gene_mp: transcript to gene id mapping <transcript:gene_name>
    @type tr_gene_mp: dict  
    @args boundary: flanking region to the signal position
    @type boundary: int
    @args sample: number of minus regions to be selected from a true region.
    @type sample: int  
    """

    #out_pos_fh = open(out_path + "/" + "tis_sig_plus_label.fa", 'w')
    out_min_fh = open("cdsstop_sig_minus_label.fa", 'w')
    true_label = 0 

    foh = helper.open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    signal_location = cdsstop_check[tr_gene_mp[fid]]

                    if loc[1] == '+': 
                        motif_seq = rec.seq[loc[2][0]:loc[2][1]-200]

                        # get index for negative signal label sequence site 
                        for ndr, stcodon in enumerate(['TAA', 'TAG', 'TGA']):
                            #search for different stop codon index 
                            idx = [xq.start() for xq in re.finditer(re.escape(stcodon), str(motif_seq).upper())]

                            # get the false labels for randomly selected region or the available ones   
                            if len(idx) > sample:
                                idx = random.sample(idx, sample)
                            
                            prev_value = dict() 
                            for nb, xp in enumerate(idx):
                                # adjusting the coordinate to the false site - mapping the relative position 
                                rloc_min = loc[2][0]+xp

                                # this has to make sure that we are removing the correct positive label 
                                if rloc_min in signal_location:
                                    continue

                                if rloc_min in prev_value:
                                    continue
                                prev_value[rloc_min] = 0
                                
                                motif_sub_seq = rec.seq[(rloc_min-boundary)+1:(rloc_min+boundary)+2]

                                # check for sanity and consensus of the fetched sequence region 
                                if not motif_sub_seq:
                                    continue
                                if not all (XJ in 'ATCG' for XJ in str(motif_sub_seq.upper())):
                                    continue
                                if not str(motif_sub_seq[boundary-1:boundary+2]).upper() in ['TAA', 'TAG', 'TGA']:
                                    continue

                                # result to fasta out
                                fseq = SeqRecord(motif_sub_seq.upper(), id='%s%s%d' % (rec.id, loc[1], rloc_min), description='-1 %s' % fid)
                                out_min_fh.write(fseq.format("fasta"))
                                true_label += 1 

                    elif loc[1] == '-': 
                        motif_seq = rec.seq[loc[2][0]+200:loc[2][1]]
                        
                        # 8713 3' UTR end, 3 nts for stop codon - 8714, 8715, 8716. 8717 end cds region
                        #print loc[0][0]-3 # 1-based coordinates 
                        #print loc[0][0] 
                        #print rec.seq[(loc[0][0]-4):(loc[0][0]-1)] # 0-based coordinates stop codon 

                        # get index for negative signal label sequence site 
                        for ndr, stcodon in enumerate(['TTA', 'TCA', 'CTA']):
                            #search for different stop codons 
                            idx = [xq.start() for xq in re.finditer(re.escape(stcodon), str(motif_seq).upper())]

                            # get the false labels for randomly selected region or the available ones   
                            if len(idx) > sample:
                                idx = random.sample(idx, sample)

                            prev_value = dict() 
                            for nb, xp in enumerate(idx):
                                # adjusting the coordinate to the false site 
                                rloc_min = loc[2][0]+xp
                                
                                #removing the true signal sequence site from selected false sites
                                if rloc_min+4 in signal_location:
                                    continue 

                                if rloc_min in prev_value:
                                    continue
                                prev_value[rloc_min] = 0 

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
                                fseq = SeqRecord(motif_sub_seq.upper(), id='%s%s%d' % (rec.id, loc[1], rloc_min), description='-1 %s' % fid)
                                out_min_fh.write(fseq.format("fasta"))
                                true_label += 1 
    out_min_fh.close()
    foh.close()
    return true_label


def false_tis_seq_fetch(fnam, Label, tis_check, tr_gene_mp, boundary=100, sample=2):
    """
    fetch the minus TIS signal label sequences 

    @args fnam: genome sequence in fasta format 
    @type fnam: str
    @args Label: signal sequence in the genome <chrom<transcript:(location, strand)>>
    @type Label: dict 
    @args tis_check: tis signals from different transcripts of the same gene <transcript:[exon_end_locations]>  
    @type tis_check: defaultdict(list)
    @args tr_gene_mp: transcript to gene id mapping <transcript:gene_name>
    @type tr_gene_mp: dict  
    @args boundary: flanking region to the signal position
    @type boundary: int
    @args sample: number of minus regions to be selected from a true region.
    @type sample: int  
    """

    out_min_fh = open("tis_sig_neg_example.fa", 'w')
    true_label = 0 

    foh = helper.open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    signal_location = tis_check[tr_gene_mp[fid]]

                    if loc[1] == '+': 
                        motif_seq = rec.seq[loc[2][0]+200:loc[2][1]]

                        # get index for negative signal label sequence site 
                        idx = [xq.start() for xq in re.finditer(re.escape('ATG'), str(motif_seq).upper())]

                        if len(idx) > sample:# limit to take maximum  false labels from one defined feature
                            idx = random.sample(idx, sample)
                        
                        prev_value = dict() 
                        # get the false labels for randomly selected region or the available ones   
                        for ndr, xp in enumerate(idx):
                            # adjusting the coordinate to the false site 
                            rloc_min = loc[2][0]+xp

                            # removing the true signl sequence site from selected false sites
                            if rloc_min+1 in signal_location:
                                continue

                            if rloc_min in prev_value:
                                continue
                            prev_value[rloc_min] = 0 

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
                            fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, loc[1], rloc_min), description='-1 %s' % fid)
                            out_min_fh.write(fseq.format("fasta"))
                            true_label += 1 

                    elif loc[1] == '-': 
                        motif_seq = rec.seq[loc[2][0]:loc[2][1]-200]

                        # get index for negative signal label sequence site 
                        idx = [xq.start() for xq in re.finditer(re.escape('CAT'), str(motif_seq).upper())]

                        if len(idx) > sample:# limit to take maximum  false labels from one defined feature
                            idx = random.sample(idx, sample)

                        prev_value = dict() 
                        # get the false labels for randomly selected region or the available ones   
                        for ndr, xp in enumerate(idx):
                            # adjusting the coordinate to the false site 
                            rloc_min = loc[2][0]+xp

                            # removing the true signal sequence site from selected false sites
                            if rloc_min+3 in signal_location:
                                continue

                            if rloc_min in prev_value:
                                continue
                            prev_value[rloc_min] = 0 
                            
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
                            fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, loc[1], rloc_min), description='-1 %s' % fid)
                            out_min_fh.write(fseq.format("fasta"))
                            true_label += 1 
    out_min_fh.close()
    foh.close()
    return true_label


def true_cdsStop_seq_fetch(fnam, Label, boundary=100):
    """
    fetch positive cdsStop signal labels 

    @args fnam: genome sequence in fasta format 
    @type fnam: str
    @args Label: signal sequence in the genome <chrom<transcript:(location, strand)>>
    @type Label: dict 
    @args boundary: flanking region to the signal position
    @type boundary: int
    """

    out_pos_fh = open("cdsstop_sig_plus_label.fa", 'w')
    true_label = 0 

    foh = helper.open_file(fnam)
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
                        fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, loc[1], int(loc[0])), description='+1 %s' % fid)
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
                        fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, loc[1], int(loc[0])), description='+1 %s' % fid)
                        out_pos_fh.write(fseq.format("fasta"))
                        true_label += 1
    out_pos_fh.close()
    foh.close()
    return true_label


def true_tis_seq_fetch(fnam, Label, boundary=100):
    """
    fetch the plus TIS signal sequence.

    @args fnam: genome sequence in fasta format 
    @type fnam: str
    @args Label: signal sequence in the genome <chrom<transcript:(location, strand)>>
    @type Label: dict 
    @args boundary: flanking region to the signal position
    @type boundary: int
    """

    out_pos_fh = open("tis_sig_pos_example.fa", 'w')
    true_label = 0 

    foh = helper.open_file(fnam)
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
                        fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, loc[1], int(loc[0])), description='+1 %s' % fid)
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
                        fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, loc[1], int(loc[0])), description='+1 %s' % fid)
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
    trans_gene_map = dict() 
    anno_db = defaultdict(list) 
    signal_point = defaultdict(list)
    
    for feature in gtf_content: # gene list from GFFParse function 
        mod_anno_db = dict()

        if signal == 'tss':
            for xp, ftid in enumerate(feature['transcripts']):
                if feature['exons'][xp].any():
                    feat_cnt += 1
                    signal_point[feature['name']].extend(feature['tss'][xp])
                    trans_gene_map[ftid[0]] = feature['name']
                    mod_anno_db[ftid[0]] = (feature['tss'][xp], 
                                feature['strand'], 
                                (int(feature['start']), int(feature['stop']))
                                )
                    break  ## considering only one transcript from a region 
        elif signal == 'tis':
            for xp, ftid in enumerate(feature['transcripts']):
                if feature['cds_exons'][xp].any():
                    feat_cnt += 1
                    signal_point[feature['name']].extend(feature['tis'][xp])
                    trans_gene_map[ftid[0]] = feature['name']
                    mod_anno_db[ftid[0]] = (feature['tis'][xp],
                                feature['strand'],
                                (int(feature['start']), int(feature['stop']))
                                )
        elif signal == 'cdsstop':
            for xp, ftid in enumerate(feature['transcripts']):
                if feature['cds_exons'][xp].any():
                    feat_cnt += 1
                    signal_point[feature['name']].extend(feature['cdsStop'][xp])
                    trans_gene_map[ftid[0]] = feature['name']
                    mod_anno_db[ftid[0]] = (feature['cdsStop'][xp],
                                feature['strand'],
                                (int(feature['start']), int(feature['stop']))
                                )
        elif signal == 'cleave':
            for xp, ftid in enumerate(feature['transcripts']):
                if feature['exons'][xp].any():
                    feat_cnt += 1
                    signal_point[feature['name']].extend(feature['cleave'][xp])
                    trans_gene_map[ftid[0]] = feature['name']
                    mod_anno_db[ftid[0]] = (feature['cleave'][xp], 
                                feature['strand'],
                                (int(feature['start']), int(feature['stop']))
                                )
        elif signal == 'splice':# considering each annotated transcript 
            for xp, ftid in enumerate(feature['transcripts']):
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
                    signal_point[feature['name']].extend(exon_ends) # appending splice sites of transcripts to the gene
                    trans_gene_map[ftid[0]] = feature['name'] # transcript to gene mapping 

        if mod_anno_db:
            anno_db[feature['chr']].append(mod_anno_db)
    
    return dict(anno_db), feat_cnt, signal_point, trans_gene_map 


def false_ss_seq_fetch(fnam, Label, don_acc_check, tr_gene_mp, boundary=100, sample=2):
    """
    false splice signals

    @args fnam: genome sequence in fasta format 
    @type fnam: str
    @args Label: signal sequence in the genome <chrom<transcript:(location, strand)>>
    @type Label: dict 
    @args don_acc_check: splice signals from different transcripts of the same gene <transcript:[exon_end_locations]>  
    @type don_acc_check: defaultdict(list)
    @args tr_gene_mp: transcript to gene id mapping <transcript:gene_name>
    @type tr_gene_mp: dict  
    @args boundary: flanking region to the signal position
    @type boundary: int
    @args sample: number of minus regions to be selected from a true region.
    @type sample: int  
    """

    don_min_fh = open("don_sig_minus_label.fa", 'w')
    acc_min_fh = open("acc_sig_minus_label.fa", 'w')

    true_label_acc = 0 
    true_label_don = 0 

    foh = helper.open_file(fnam)
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

                            if len(idx) > sample: # limit to take maximum 2 false labels from one defined feature
                                idx = random.sample(idx, sample)
                            
                            prev_value = dict()
                            for ndr, xp in enumerate(idx):
                                rloc_min = int(signal_location[0])+xp # adjusting the coordinate to the false site 
                                if rloc_min+3 in signal_location: # ---AG|Exon~~~~ removing the true signal from false site 
                                    continue

                                if rloc_min in prev_value:
                                    continue

                                prev_value[rloc_min]=0 

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
                                fseq_acc = SeqRecord(acc_mot_seq.upper(), id='%s%s%d' % (rec.id, loc[-1], rloc_min), description='-1 %s' % fid)
                                acc_min_fh.write(fseq_acc.format("fasta"))
                                true_label_acc += 1 
                        
                        if not numpy.isnan(loc[1]):
                            # donor splice site signal 
                            don_t_seq = rec.seq[int(signal_location[0]):int(signal_location[-1])] # from gene body 
                            # get index for negative signal label sequence site 
                            idx = [xq.start() for xq in re.finditer(re.escape('GT'), str(don_t_seq).upper())]

                            if len(idx) > sample:# limit to take maximum 2 false labels from one defined feature
                                idx = random.sample(idx, sample)

                            prev_value = dict()
                            for ndr, xj in enumerate(idx):
                                rloc_pos = int(signal_location[0])+xj# adjusting the coordinate to the false site 
                                if rloc_pos in signal_location:
                                    continue

                                if rloc_pos in prev_value:
                                    continue

                                prev_value[rloc_pos]=0 

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
                                fseq_don = SeqRecord(don_mot_seq.upper(), id='%s%s%d' % (rec.id, loc[-1], rloc_pos), description='-1 %s' % fid)
                                don_min_fh.write(fseq_don.format("fasta"))
                                true_label_don += 1 

                    elif loc[-1]=="-":
                        if not numpy.isnan(loc[0]):
                            # donor splice signal site 
                            don_t_seq = rec.seq[int(signal_location[0]):int(signal_location[-1])]
                            # get index for negative signal label sequence site 
                            idx = [xq.start() for xq in re.finditer(re.escape('AC'), str(don_t_seq).upper())]

                            if len(idx) > sample:# limit to take maximum 2 false labels from one defined feature
                                idx = random.sample(idx, sample)

                            prev_value = dict()
                            for ndr, xp in enumerate(idx):
                                # adjusting the coordinate to the false site 
                                rloc_min = int(signal_location[0])+xp
                                if rloc_min+3 in signal_location:
                                    continue

                                if rloc_min in prev_value:
                                    continue

                                prev_value[rloc_min] = 0 

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
                                fseq_don = SeqRecord(don_mot_seq.upper(), id='%s%s%d' % (rec.id, loc[-1], rloc_min), description='-1 %s' % fid)
                                don_min_fh.write(fseq_don.format("fasta"))
                                true_label_don += 1 
                        
                        if not numpy.isnan(loc[1]):
                            # acceptor splice signal site 
                            acc_t_seq = rec.seq[int(signal_location[0]):int(signal_location[-1])]
                            # get index for negative signal label sequence site 
                            idx = [xq.start() for xq in re.finditer(re.escape('CT'), str(acc_t_seq).upper())]

                            if len(idx) > sample:# limit to take maximum 2 false labels from one defined feature
                                idx = random.sample(idx, sample)

                            prev_value = dict()
                            for ndr, xk in enumerate(idx):
                                # adjusting the coordinate to the false site 
                                rloc_pos = int(signal_location[0])+xk
                                if rloc_pos in signal_location:
                                    continue

                                if rloc_pos in prev_value:
                                    continue
                                
                                prev_value[rloc_pos] = 0 

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
                                fseq_acc = SeqRecord(acc_mot_seq.upper(), id='%s%s%d' % (rec.id, loc[-1], rloc_pos), description='-1 %s' % fid)
                                acc_min_fh.write(fseq_acc.format("fasta"))
                                true_label_acc += 1

    foh.close()
    don_min_fh.close()
    acc_min_fh.close()
    return true_label_acc, true_label_don


def minus_tss_seq_fetch(fnam, Label, tss_check, tr_gene_mp, boundary=100, sample=3):
    """
    fetch the minus TSS signal sequence label

    @args fnam: genome sequence in fasta format 
    @type fnam: str
    @args Label: signal sequence in the genome <chrom<transcript:(location, strand)>>
    @type Label: dict 
    @args tss_check: tss signals from different transcripts of the same gene <transcript:[exon_end_locations]>  
    @type tss_check: defaultdict(list)
    @args tr_gene_mp: transcript to gene id mapping <transcript:gene_name>
    @type tr_gene_mp: dict  
    @args boundary: flanking region to the signal position
    @type boundary: int
    @args sample: number of minus regions to be selected from a true region.
    @type sample: int  
    """

    out_min_fh = open("tss_sig_neg_example.fa", 'w')
    true_label = 0 

    foh = helper.open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    signal_location = tss_check[tr_gene_mp[fid]]

                    prev_value = dict() 
                    for ndr in range(sample):

                        if loc[1]=='+': # I---------|~~~~~
                            rloc = random.randint(loc[2][0]-1000,loc[2][1])
                        elif loc[1]=='-':
                            rloc = random.randint(loc[2][0],loc[2][1]+1000)
                        else: ## if there is no strand info 
                            rloc = random.randint(loc[2][0]-1000,loc[2][1])

                        # remove the true signal index from random sampling 
                        if rloc in signal_location:
                            continue

                        if rloc in prev_value:
                            continue
                        prev_value[rloc] = 0  

                        #motif_seq = rec.seq[rloc-boundary:rloc+boundary+1]
                        motif_seq = rec.seq[rloc-boundary:rloc+boundary]

                        # sanity check for featched sequence 
                        if len(motif_seq) != 2*boundary:
                            continue
                        if not motif_seq:
                            continue
                        if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                            continue

                        # write to fasta out 
                        fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, loc[1], rloc), description='-1 %s' % fid)
                        out_min_fh.write(fseq.format("fasta"))
                        true_label += 1 
    out_min_fh.close()
    foh.close()
    return true_label


def minus_cleave_seq_fetch(fnam, Label, cleave_check, tr_gene_mp, boundary=100, sample=2):
    """
    fetch the minus TSS signal sequence label

    @args fnam: genome sequence in fasta format 
    @type fnam: str
    @args Label: signal sequence in the genome <chrom<transcript:(location, strand)>>
    @type Label: dict 
    @args cleave_check: tss signals from different transcripts of the same gene <transcript:[exon_end_locations]>  
    @type cleave_check: defaultdict(list)
    @args tr_gene_mp: transcript to gene id mapping <transcript:gene_name>
    @type tr_gene_mp: dict  
    @args boundary: flanking region to the signal position
    @type boundary: int
    @args sample: number of minus regions to be selected from a true region.
    @type sample: int  
    """

    out_min_fh = open("cleave_sig_minus_label.fa", 'w')
    true_label = 0 

    foh = helper.open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():
                    if loc[2][1]-loc[2][0] < 200: # remove small transcripts 
                        continue

                    signal_location = cleave_check[tr_gene_mp[fid]]

                    prev_value = dict() 
                    for ndr in range(sample):
                        if loc[1]=='+': # ---------|~~~~~
                            rloc = random.randint(loc[2][0],loc[2][1]-200)
                        elif loc[1]=='-':
                            rloc = random.randint(loc[2][0]+200,loc[2][1])

                        # remove the true signal index from random sampling 
                        if rloc in signal_location:
                            continue

                        if rloc in prev_value:
                            continue
                        prev_value[rloc] = 0 

                        #motif_seq = rec.seq[rloc-boundary:rloc+boundary+1]
                        motif_seq = rec.seq[rloc-boundary:rloc+boundary]

                        # sanity check for featched sequence 
                        #if len(motif_seq) != 2*boundary+1:
                        if len(motif_seq) != 2*boundary:
                            continue
                        if not motif_seq:
                            continue
                        if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                            continue

                        # write to fasta out 
                        fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, loc[1], rloc), description='-1 %s' % fid)
                        out_min_fh.write(fseq.format("fasta"))
                        true_label += 1 
    out_min_fh.close()
    foh.close()
    return true_label


def true_ss_seq_fetch(fnam, Label, boundary=100):
    """
    True splice signals 

    @args fnam: genome sequence in fasta format 
    @type fnam: str
    @args Label: random signal label position in the genome  <chrom<transcript:(location, strand)>>
    @type Label: dict
    @args boundary: flanking region to the signal upstream and downstream
    @type boundary: integer 
    """

    don_pos_fh = open("don_sig_plus_label.fa", 'w')
    acc_pos_fh = open("acc_sig_plus_label.fa", 'w')

    true_label_acc = 0
    true_label_don = 0 

    foh = helper.open_file(fnam)
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
                            fseq_acc = SeqRecord(acc_mot_seq.upper(), id='%s%s%d' % (rec.id, loc[-1], int(loc[0])), description='+1 %s' % fid)
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
                            fseq_don = SeqRecord(don_mot_seq.upper(), id='%s%s%d' % (rec.id, loc[-1], int(loc[1])), description='+1 %s' % fid)
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
                            fseq_don = SeqRecord(don_mot_seq.upper(), id='%s%s%d' % (rec.id, loc[-1], int(loc[0])), description='+1 %s' % fid)
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
                            fseq_acc = SeqRecord(acc_mot_seq.upper(), id='%s%s%d' % (rec.id, loc[-1], int(loc[1])), description='+1 %s' % fid)
                            acc_pos_fh.write(fseq_acc.format("fasta"))
                            true_label_acc += 1 
    don_pos_fh.close()
    acc_pos_fh.close()
    foh.close()
    return true_label_acc, true_label_don


def plus_tss_cleave_seq_fetch(signal, fnam, Label, boundary=100):
    """
    fetch the plus TSS and Cleave signal sequence

    @args signal: signal type cleave or tss 
    @type signal: str
    @args fnam: genome sequence in fasta format 
    @type fnam: str
    @args Label: signal sequence in the genome <chrom<transcript:(location, strand)>>
    @type Label: dict 
    @args boundary: flanking region to the signal position
    @type boundary: int
    """

    out_pos_fh = open(signal + "_sig_pos_example.fa", 'w')
    true_label = 0 

    foh = helper.open_file(fnam)
    for rec in SeqIO.parse(foh, "fasta"):
        if rec.id in Label:
            for Lsub_feat in Label[rec.id]:
                for fid, loc in Lsub_feat.items():

                    #motif_seq = rec.seq[int(loc[0])-boundary:int(loc[0])+boundary+1] # TSS ---A---
                    motif_seq = rec.seq[int(loc[0])-boundary:int(loc[0])+boundary] # TSS ---A--

                    if loc[1] == '-': 
                        motif_seq = motif_seq.reverse_complement()

                    # sanity check for the fetched sequence 
                    if len(motif_seq) != boundary*2: 
                        continue
                    if not motif_seq:
                        continue
                    if not all (XJ in 'ATCG' for XJ in str(motif_seq.upper())):
                        continue

                    # write to fasta out 
                    fseq = SeqRecord(motif_seq.upper(), id='%s%s%d' % (rec.id, loc[1], int(loc[0])), description='+1 %s' % fid)
                    out_pos_fh.write(fseq.format("fasta"))
                    true_label += 1 

    out_pos_fh.close()
    foh.close()
    return true_label


def select_labels(feat_db, feat_count, label_cnt):
    """
    Random sampling to select signal lables

    @args feat_db: signal sequence location from the genome 
    @type feat_db: dict
    @args feat_count: annotated total signal features 
    @type feat_count: integer
    @args label_cnt: number of features to be considered
    @type label_cnt: integer
    """

    assert label_cnt <= feat_count, 'Number of features annotated %d' % feat_count

    try:
        accept_prob = (1.0*label_cnt)/feat_count
    except:
        accept_prob = 1
    
    while True: # ensure the label count 
        counter, LSet = recursive_fn(feat_db, label_cnt, accept_prob)
        if label_cnt <= counter:
            break
        sys.stdout.write('    still trying ... %d\n' % counter) 

    return LSet, counter


def recursive_fn(f_db, lb_cnt, apt_prob):
    """
    This function returns the random samples based on the label counts 

    @args f_db: signal sequence location from genome 
    @type f_db: dict 
    @args lb_cnt: number of labels 
    @type lb_cnt: integer
    @args apt_prob: accept probability to make a the random search  entire list
    @type apt_prob: float 
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
        fas_name = sys.argv[1] # genome sequence - fasta 
        gff_name = sys.argv[2] # genome annotation - gff/gtf 
    except:
        print __doc__
        sys.exit(-1) 
    
    get_feature_seq(fas_name, gff_name)
