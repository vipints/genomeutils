#!/usr/bin/env python 
"""
modules for large scale experiment runs, example: generating standard genome indices. 

Usage: 
import prepare_data as pd 
pd.create_star_genome_index.__doc__

Requirement:
    STAR aligner: https://github.com/alexdobin/STAR 
    gfftools: https://github.com/vipints/genomeutils/tree/master/gfftools 
"""

import os 
import re
import sys 
import shutil
import subprocess 
from Bio import SeqIO 
from gfftools import helper, GFFParser


def stop_err(msg):
    """
    stop the execution and print out the captured error message. 
    """
    sys.stderr.write('%s\n' % msg)
    sys.exit(-1)


def clean_genome_file(chr_names, fas_file, fas_out):
    """
    make a stable genome file with valid contigs 

    @args chr_names: different contig names with a valid genome sequence 
    @type chr_names: dict
    @args fas_file: genome sequence in fasta file  
    @type fas_file: str 
    @args fas_out: new genome sequence file in fasta format 
    @type fas_out: str 
    """

    # get the valid chromosome identifier from user as STDIN
    chr_names = dict()
    
    # get the filehandler from input file
    try:
        fh = helper.open_file(fas_file)
    except Exception, errmsg:
        stop_err('error in reading file '+ errmsg) 

    # check the out filehandler
    try:
        outfh = open(fas_out, "w")
    except Exception, errmsg:
        stop_err('error in writing file '+ errmsg) 

    # writing stable contig genome sequence in FASTA format 
    for rec in SeqIO.parse(fh, "fasta"):
        if rec.id in chr_names:
            outfh.write(rec.format("fasta"))

    fh.close()
    outfh.close()


def make_anno_db(gff_file): 
    """
    extract the features from a gtf/gff file and store efficiently to query 

    @args gff_file: genome annotation file
    @type gff_file: str 
    """

    gff_cont = GFFParser.Parse(gff_file)  

    intron_size = dict() 
    exon_size = dict() 

    for rec in gff_cont:
        for idx, tid in enumerate(rec['transcripts']):

            if not rec['exons'][idx].any():
                continue

            try: # (Pdb) rec['exons'][0] -> array(nan)
                import numpy as np 
                if np.isnan(rec['exons'][idx]):
                    continue
            except:
                pass 
                    
            try:
                exon_cnt = len(rec['exons'][idx])
            except:
                import pdb 
                pdb.set_trace()

            if exon_cnt > 1:
                intron_start = 0 
                
                for xq, excod in enumerate(rec['exons'][idx]): 
                    
                    if xq > 0: 
                        #print intron_start, excod[0]-1 
                        if excod[0]-intron_start==1:
                            intron_start = excod[1]+1
                            exon_size[intron_start-excod[0]] = 1
                            continue

                        intron_size[excod[0]-intron_start] = 1 
                        #print excod[0]-intron_start

                    intron_start = excod[1]+1
                    exon_size[intron_start-excod[0]] = 1

                    #print intron_start-excod[0]
    if intron_size:
        keys_int = sorted(intron_size)
        print 'MinIntronLength %d %d %d'  %(keys_int[0], keys_int[1], keys_int[2])
        print 'MaxIntronLength %d %d %d'  %(keys_int[-1], keys_int[-2], keys_int[-3])
        print 
        keys_ex = sorted(exon_size)
        print 'MinExonLength %d %d %d'  %(keys_ex[0], keys_ex[1], keys_ex[2]) 
        print 'MaxExonLength %d %d %d'  %(keys_ex[-1], keys_ex[-2], keys_ex[-3]) 
    else:
        print "Error in feature mapping, please check the source of parent child features" 
        print "May be the sources are different for parents and child features of the parent Gene"


def create_star_genome_index(fasta_file, out_dir, genome_anno=None, num_workers=1, onematelength=100):
    """
    Creating STAR genome index with or without using genome annotation

    @args fasta_file: reference genome sequence file .fasta format 
    @type fasta_file: str 
    @args out_dir: genome index binary file storage place  
    @type out_dir: str 
    @args genome_anno: genome annotation file (optional) 
    @type genome_anno: str 
    @args num_workers: number of threads to run (default value = 1)
    @type num_workers: int 
    @args onematelength: One Mate Length (default value=100) 
    @type num_workers: int 
    """
    
    if not genome_anno:
        cli_cmd = 'STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --runThreadN %d' % (out_dir, fasta_file, num_workers) 
    else:
        ## check for the file type  
        gff_hand = helper.open_file(genome_anno)
    
        for rec in gff_hand:
            rec = rec.strip('\n\r')

            # skip empty line fasta identifier and commented line
            if not rec or rec[0] in  ['#', '>']:
                continue
            # skip the genome sequence 
            if not re.search('\t', rec):
                continue

            parts = rec.split('\t')
            assert len(parts) >= 8, rec

            ftype, tags = GFFParser.attribute_tags(parts[-1])
            break 

        gff_hand.close() 

        ## according to the file type 
        if ftype:
            cli_cmd = 'STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --runThreadN %d --sjdbGTFfile %s --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang %d' % (out_dir, fasta_file, num_workers, genome_anno, onematelength) 
        else:
            cli_cmd = 'STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --runThreadN %d --sjdbGTFfile %s --sjdbOverhang %d' % (out_dir, fasta_file, num_workers, genome_anno, onematelength) 

    ## create downloadpath if doesnot exists 
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    ## start the index job 
    try:
        process = subprocess.Popen(cli_cmd, shell=True) 
        process.wait()
    except:
        print "error"
        sys.exit(-1)

    print 
    print "STAR genome index files are stored at %s" % out_dir
    print 


if __name__=="__main__":
    print __doc__

"""
    fasta_file = "/home/data/C_elegans/ce10.fa"
    out_dir = "/home/data/genome/"
    genome_anno = "/home/data/ce10.gff"

    create_star_genome_index(fasta_file, out_dir, genome_anno, num_workers=1, onematelength=100):
"""
