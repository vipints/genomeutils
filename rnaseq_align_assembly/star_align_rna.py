#!/usr/bin/env python 
"""
Program to align rnaseq transcriptome reads to the genome using STAR aligner 

Requirement: 
    STAR - https://github.com/alexdobin/STAR/releases
    pysam - http://code.google.com/p/pysam/ 
"""

import os 
import sys 
import subprocess

def run_star_alignment(org_db, read_type='PE', max_mates_gap_length=10000, num_cpus=1):
    """
    wrapper for running STAR program 

    @args org_db: a python dictionary with all details about a single organism 
    @type org_db: dict

    """

    #os.environ['PATH'] += os.pathsep + '/home/share/software/STAR_2.3.0e.Linux_x86_64/'

    genome_dir = org_db['index']
    gtf_db = org_db['gtf']

    if read_type == 'PE':
        read_file = "%s %s" % (org_db['fastq'][0], org_db['fastq'][1])
    else:
        read_file = org_db['fastq'][0]
    
    ## getting the command to uncompress the read file, store the files in compressed form  
    zip_type = {".gz" : "zcat", ".bz2" : "bzcat"} 
    file_prefx, ext = os.path.splitext(org_db['fastq'][0])

    max_lenth_intron = org_db['max_intron']

    make_file_name = "mkfifo Aligned.out.sam"
    process = subprocess.Popen(make_file_name, shell=True) 
    make_bg_process = "cat Aligned.out.sam | samtools view -Shb - | samtools sort - Aligned.sort.out &"
    process = subprocess.Popen(make_bg_process, shell=True) 

    make_star_run = "STAR \
    --genomeDir %s \
    --readFilesIn %s \
    --readFilesCommand %s \
    --runThreadN %d \
    --outFilterMultimapScoreRange 2 \
    --outFilterMultimapNmax 100 \
    --outFilterMismatchNmax 10 \
    --alignIntronMax %d \
    --sjdbGTFfile %s \
    --sjdbScore 1 \
    --sjdbOverhang 5 \
    --genomeLoad LoadAndRemove" % (genome_dir, read_file, zip_type[ext], num_cpus, max_lenth_intron, gtf_db)

    sys.stdout.write('\trunning STAR program as: %s \n' % make_star_run)
    process = subprocess.Popen(make_star_run, shell=True) 
    process.wait()


def uniq_mapped_reads(bam_file, multi_map=1):
    """
    Filter a BAM file to extract only the uniquely mapped reads

    @args bam_file: binary file for storing the sequencing reads.
    @type bam_file: str 
    @args multi_map: number of hits of a read (default 1) 
    @type multi_map: integer 
    """

    import pysam

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


def calculate_insert_size_from_bam(bam_file):
    """
    estimate the insert size from uniquely mapped reads in a BAM file 
    """

    import pysam 


def calculate_insert_size(org_db):
    """
    calculate the insert-size from raw read sequence file
    """

    import re 
    from glob import glob 


    #FIXME 
    """
    vipin@gpu-3-9: /cbio/grlab/nobackup/SignalPrediction/SRA-rnaseq/S_enterica/source_data$ python ~/tmp/7281991/estimate-insert-sizes /cbio/grlab/share/databases/genomes/S_enterica/ensembl_release-21/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.GCA_000006945.1.21.dna.toplevel.fa SRR863221_1.fastq.bz2 SRR863221_2.fastq.bz2
    Processing:
     SRR863221_1.fastq.bz2
      SRR863221_2.fastq.bz2
      [M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (0, 0, 0, 0)
      [M::mem_pestat] skip orientation FF as there are not enough pairs
      [M::mem_pestat] skip orientation FR as there are not enough pairs
      [M::mem_pestat] skip orientation RF as there are not enough pairs
      [M::mem_pestat] skip orientation RR as there are not enough pairs
      Traceback (most recent call last):
        File "/cbio/grlab/home/vipin/tmp/7281991/estimate-insert-sizes", line 100, in <module>
            mean = most_likely[1]['mean']
            KeyError: 'mean'

    """
