#!/usr/bin/env python 
"""
master script to execute the pipeline
"""

import os 
import sys 
import yaml 

import libpyjobrunner as pg

from optparse import OptionParser

from fetch_remote_data import download_data as dld
from fetch_remote_data import prepare_data as ppd

from signal_labels import experiment_details_db as expdb

assert sys.version_info[:2] >= ( 2, 4 )

def main():
    """
    Managing the experiment run in different levels 

    Options
    
    -1 download_public_data from SRA ENSEMBL Phytozome
        uncompressing the SRA file 
        manual cleaning of downloaded genome 
        create the genome genome indices

    -2 rnaseq read mapping process     
    """

    parser = OptionParser() 

    parser.add_option( "-1", "--download_public_data", action="store_true", dest="download_public_data", default=False, help="download public datasets")
    parser.add_option( "-a", "--genome_index", action="store_true", dest="genome_index", default=False, help="Create STAR genome index to align the reads.")
    parser.add_option( "-2", "--read_mapping", action="store_true", dest="read_mapping", default=False, help="RNASeq read mapping to the genome using STAR.")

    ( options, args ) = parser.parse_args()
    try:
        config_file = args[0]
    except:
        print __doc__
        sys.exit(-1)

    print 'Using config file %s for the experiment.' % config_file

    if not (options.download_public_data ^ options.genome_index ^ \
            options.read_mapping):
        parser.print_help()
        sys.exit(-1)

    if options.download_public_data:
        print 'Operation selected: Download public genome dataset and Sequencing experiment files'
        download_public_data(config_file)

    if options.genome_index:
        print 'Operation selected: Create STAR genome index'
        create_genome_index(config_file)

    if options.read_mapping: 
        print 'Operation selected: Read alignment with STAR'
        align_rnaseq_reads(config_file)


def call_align_reads(args_list):
    """
    wrapper for submitting jobs to pygrid
    """

    from rnaseq_align_assembly import star_align_rna as rnastar 
    org_db, read_type, max_mates_gap_length, num_cpus = args_list
    rnastar.run_star_alignment(org_db, read_type, max_mates_gap_length, num_cpus) 
    return 'done'


def align_rnaseq_reads(yaml_config):
    """
    wrapper for aligning rnaseq reads using 
    """

    orgdb = expdb.experiment_db(yaml_config)
    Jobs = []
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        lib_type = 'PE'
        lib_type = 'SE' if len(det['fastq'])==1 else lib_type

        arg = [[det, lib_type, 100000, 4]]

        job = pg.cBioJob(call_align_reads, arg) 
    
        job.mem="32gb"
        job.vmem="32gb"
        job.pmem="8gb"
        job.pvmem="8gb"
        job.nodes = 1
        job.ppn = 4
        job.walltime = "08:00:00"
        
        Jobs.append(job)

    print 
    print "sending jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


def call_genome_index(args_list):
    """
    wrapper for submitting jobs to pygrid
    """

    fasta_file, out_dir, genome_anno, num_workers, onematelength = args_list
    ppd.create_star_genome_index(fasta_file, out_dir, genome_anno, num_workers, onematelength)
    return 'done'


def create_genome_index(yaml_config):
    """
    wrapper for calling genome index function 
    """

    orgdb = expdb.experiment_db(yaml_config)

    Jobs = []
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        arg = [[det['fasta'], det['genome_index_dir'], det['gtf'], 1, det['read_length']-1]]

        job = pg.cBioJob(call_genome_index, arg) 
    
        job.mem="8gb"
        job.vmem="8gb"
        job.pmem="8gb"
        job.pvmem="8gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "01:00:00"
        
        Jobs.append(job)

    print 
    print "sending jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


def download_public_data(yaml_config):
    """
    """

    config_map = yaml.safe_load(open(yaml_config, "rU"))
    #import ipdb 
    #ipdb.set_trace()

    print config_map


def download_fasta(org_details):
    """
    download fasta file from remote data publishing services
    """

    for org_name, det in org_details.items():
        if det['release_db'] == 'ensembl':
            dld.fetch_ensembl_fasta(det['release_num'], det['name'], det['fasta'])
        elif det['release_db'] == 'phytozome':
            dld.fetch_phytozome_fasta(det['release_num'], det['name'], det['fasta']) 
        elif det['release_db'] == 'ensembl_metazoa':
            dld.fetch_ensembl_metazoa_fasta(det['release_num'], det['name'], det['fasta'])
        else:
            print "download fasta plugin for %s not available, module works with ensembl, ensembl_metazoa and phytozome." % det['release_db']


def download_gtf(org_details):
    """
    download gtf/gff file from remote data publishing services
    """
    for org_name, det in org_details.items():
        if det['release_db'] == 'ensembl':
            dld.fetch_ensembl_gtf(det['release_num'], det['name'], det['gtf'])
        elif det['release_db'] == 'phytozome':
            dld.fetch_phytozome_gff(det['release_num'], det['name'], det['gtf']) 
        elif det['release_db'] == 'ensembl_metazoa':
            dld.fetch_ensembl_metazoa_gtf(det['release_num'], det['name'], det['gtf'])
        else:
            print "download gtf plugin for %s not available, module works with ensembl, ensembl_metazoa and phytozome." % det['release_db']


def download_uncompress_sra_file(org_details):
    """
    download and uncompress the sra files 
    """

    for org_name, det in org_details.items():
        sra_file = dld.download_sra_file(det['sra_run_id'], det['fastq_path'])
        
        library_type = "pe"
        compress_format = "gzip"
        dld.uncompress_sra_file(sra_file, det['fastq_path'], library_type, compress_format)


if __name__=="__main__":
    
    main() 

    """
    from signal_labels import org_details_db as odb
    org_details = odb.make_org_db(infile, data_path, exp_path) 

    download_fasta(org_details) 

    download_gtf(org_details)     

    download_uncompress_sra_file(org_details)

    #TODO
    save a pickle file with organisms details with updated genome annotation and sra file, path informations 

    the object will be passed to the next preprocessing manual tweeking

    
    #FIXME 
    fas_file = ""
    chr_names = prd.read_genome_file(fas_file) 

    fas_out = "" 
    prd.clean_genome_file(chr_names, fas_file, fas_out):

    gtf_out = "" 
    prd.clean_anno_file(chr_names, gtf_file, gtf_out)
    """
