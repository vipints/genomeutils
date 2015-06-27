#!/usr/bin/env python 
"""
Program to align rnaseq transcriptome reads to the genome using STAR aligner 

Requirement: 
    STAR - https://github.com/alexdobin/STAR/releases
    pysam - http://code.google.com/p/pysam/ 
    mmr - https://github.com/ratschlab/mmr 
"""

import os 
import re
import sys 
import subprocess
from collections import defaultdict


def run_mmr(org_name, read_map_dir, threads=3):
    """
    a pythonic wrapper for multiple mapper resolution program

    @args org_name: Organism name, example case A_thaliana 
    @type org_name: str 
    @args read_map_dir: directory where the STAR bam (aligned reads) file located
    @type read_map_dir: str 
    @args threads: number of threads to use for the run (default: 3)
    @type threads: int  
    """
    import pysam

    try:
        subprocess.call(["mmr"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        exit("Please make sure that the `mmr` binary is in your $PATH")

    ## mmr works well with bam file sorted by read id  
    bam_file = "%s/%s_Aligned.sortedByName.out.bam" % (read_map_dir, org_name) 
    if not os.path.isfile(bam_file):
        sys.stdout.write("warning: failed to fetch read id sorted BAM file for organism: %s, trying to get the raw alignment file\n" % org_name)

        bam_file = "%s/%s_Aligned.out.bam" % (read_map_dir, org_name) ## unsorted bam file from STAR output  
        if not os.path.isfile(bam_file):
            exit("error: failed to fetch STAR read alignment file for %s %s\n" % (org_name, bam_file))

        ## sorting bam file 
        sorted_bam = "%s/%s_Aligned.sortedByName.out" % (read_map_dir, org_name)
        if not os.path.isfile("%s.bam" % sorted_bam):
            sys.stdout.write("trying to sort based by read id with output prefix as: %s\n" % sorted_bam)
            pysam.sort("-n", bam_file, sorted_bam)

        bam_file = "%s.bam" % sorted_bam
    
    sys.stdout.write("using bam file from %s\n" % bam_file)
    outFile = "%s/%s_Aligned_mmr.bam" % (read_map_dir, org_name) 

    iterations = 3 
    ## provide a bam file sorted by read id
    cli_mmr = "module load gcc; mmr -b -p -V -t %d -I %d -o %s %s" % (threads, iterations, outFile, bam_file)  

    try:
        sys.stdout.write('\trun MMR as: %s \n' % cli_mmr)
        ## changing the working dir to run mmr 
        os.chdir(read_map_dir)

        process = subprocess.Popen(cli_mmr, shell=True) 
        returncode = process.wait()

        if returncode !=0:
            raise Exception, "Exit status return code = %i" % returncode

        sys.stdout.write('MMR run finished. result file stored at %s\n' % outFile)
    except Exception, e:
        exit('Error running MMR.\n%s' %  str( e ))


def run_star_alignment(org_db, read_type='PE', max_mates_gap_length=100000, num_cpus=1):
    """
    wrapper for running STAR program 

    @args org_db: a python dictionary with all details about a single organism 
    @type org_db: defaultdict
    @args read_type: library type - paired-end or single-end (default: PE)
    @type read_type: str 
    @args max_mates_gap_length: maximum insert size from the sample (default: 10000)
    @type max_mates_gap_length: int 
    @args num_cpus: number of threads to use for the run (default: 1)
    @type num_cpus: int 
    """
    try:
        subprocess.call(["STAR"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        exit("Please make sure that the `STAR` binary is in your $PATH")

    from gfftools import helper, GFFParser

    genome_dir = org_db['genome_index_dir']## genome indices and annotation file
    gtf_db = org_db['gtf']

    ## check for the annotation file type gff or gtf 
    gff_hand = helper.open_file(gtf_db)
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

    ## library type 
    if read_type == 'PE':
        read_file = "%s/%s %s/%s" % (org_db['fastq_path'], org_db['fastq'][0], org_db['fastq_path'], org_db['fastq'][1])
    else:
        read_file = "%s/%s" % (org_db['fastq_path'], org_db['fastq'][0])
    
    ## getting the command to uncompress the read file
    zip_type = {".gz" : "zcat", ".bz2" : "bzcat"} 
    file_prefx, ext = os.path.splitext(org_db['fastq'][0])

    ## genomic feature information 
    max_lenth_intron = org_db['max_intron_len']

    out_prefix = '%s/%s_' % (org_db['read_map_dir'], org_db['short_name'])

    ## according to the file type 
    if ftype:
        make_star_run = "STAR \
        --genomeDir %s \
        --readFilesIn %s \
        --readFilesCommand %s \
        --outFileNamePrefix %s \
        --runThreadN %d \
        --outFilterMultimapScoreRange 2 \
        --outFilterMultimapNmax 30 \
        --outFilterMismatchNmax 4 \
        --alignIntronMax %d \
        --sjdbGTFfile %s \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbScore 1 \
        --sjdbOverhang 5 \
        --alignMatesGapMax %d\
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM Unsorted \
        --genomeLoad LoadAndRemove" % (genome_dir, read_file, 
            zip_type[ext], out_prefix, num_cpus, max_lenth_intron, gtf_db, max_mates_gap_length)
    else:
        make_star_run = "STAR \
        --genomeDir %s \
        --readFilesIn %s \
        --readFilesCommand %s \
        --outFileNamePrefix %s \
        --runThreadN %d \
        --outFilterMultimapScoreRange 2 \
        --outFilterMultimapNmax 30 \
        --outFilterMismatchNmax 4 \
        --alignIntronMax %d \
        --sjdbGTFfile %s \
        --sjdbGTFfeatureExon exon \
        --sjdbScore 1 \
        --sjdbOverhang 5 \
        --alignMatesGapMax %d\
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM Unsorted \
        --genomeLoad LoadAndRemove" % (genome_dir, read_file, 
            zip_type[ext], out_prefix, num_cpus, max_lenth_intron, gtf_db, max_mates_gap_length)

    sys.stdout.write('\trunning STAR program as: %s \n' % make_star_run)
    try:
        process = subprocess.Popen(make_star_run, shell=True) 
        process.wait()

        print "star run completed. result file stored at %sAligned.out.bam" % out_prefix
    except Exception, e:
        print 'Error running STAR.\n%s' %  str( e )
        sys.exit(0)


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

    import time 
    t_start = time.time() 

    sam_file = pysam.Samfile(bam_file, "rb") 

    ## defining the out bam file 
    file_pref, ext = os.path.splitext(bam_file)
    uniq_map_bam_file = '%s_uniq_map_reads.bam' % file_pref
    
    try: 
        out_bam_fh = pysam.Samfile(uniq_map_bam_file, "wb", template=sam_file)
    except:
        exit("error: cannot create the file %s." % uniq_map_bam_file) 

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

    time_taken = time.time() - t_start
    print "time taken in seconds ", time_taken


def read_directions_count(bam_file):
    """
    get the reads directions count from a bam file 

    @args bam_file: binary file formt for storing sequencing reads information
    @type bam_file: str 
    """

    import pysam 

    ## indexing the in bam file 
    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file) 

    reverse_cnt = 0 
    forward_cnt = 0

    bam_fh = pysam.Samfile(bam_file, "rb") 

    for read in bam_fh.fetch():
         if read.is_proper_pair and read.is_read1:
            if read.is_reverse:
                reverse_cnt += 1
            else:
                forward_cnt += 1 

    bam_fh.close() 
    return {'forward_reads_count': forward_cnt, 'reverse_reads_count': reverse_cnt} 


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
    $ python ~/tmp/7281991/estimate-insert-sizes Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.GCA_000006945.1.21.dna.toplevel.fa SRR863221_1.fastq.bz2 SRR863221_2.fastq.bz2
    Processing:
     SRR863221_1.fastq.bz2
      SRR863221_2.fastq.bz2
      [M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (0, 0, 0, 0)
      [M::mem_pestat] skip orientation FF as there are not enough pairs
      [M::mem_pestat] skip orientation FR as there are not enough pairs
      [M::mem_pestat] skip orientation RF as there are not enough pairs
      [M::mem_pestat] skip orientation RR as there are not enough pairs
      Traceback (most recent call last):
        File "/home/vipin/tmp/7281991/estimate-insert-sizes", line 100, in <module>
            mean = most_likely[1]['mean']
            KeyError: 'mean'

    """
