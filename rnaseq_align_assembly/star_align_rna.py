#!/usr/bin/env python 
"""
Program to align rnaseq transcriptome reads to the genome using STAR aligner 

Requirement: 
    STAR - https://github.com/alexdobin/STAR/releases
    pysam - http://code.google.com/p/pysam/ 
"""

import os 
import pysam


def run_star_alignment():
    """
    """

    genome_dir = ""
    

def uniq_mapped_reads(bam_file, multi_map=1):
    """
    Filter a BAM file to extract only the uniquely mapped reads

    @args bam_file: binary file for storing the sequencing reads.
    @type bam_file: str 
    @args multi_map: number of hits of a read  
    @type multi_map: integer 
    """

    ## indexing the in bam file 
    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file) 

    #import time 
    #t_start = time.time() 

    sam_file = pysam.Samfile(bam_file, "rb") 

    ## getting the abs path of in bam file  
    bam_file_path = os.path.abspath(bam_file)
    bam_file_dir = os.path.dirname(bam_file_path)
    
    ## defining the out bam file 
    uniq_map_bam_file = '%s/unique_map_reads.bam' % bam_file_dir 
    out_bam_fh = pysam.Samfile(uniq_map_bam_file, "wb", template=sam_file)

    ## filtering the alignment based on the read hits
    for read in sam_file.fetch():
        for tag in read.tags:
            if tag[0] == 'NH':
                if int(tag[1]) <= multi_map:
                    out_bam_fh.write(read) 

    sam_file.close() 
    out_bam_fh.close() 

    ## indexing the out bam file 
    pysam.index(uniq_map_bam_file)

    #time_taken = time.time() - t_start
    #print time_taken
