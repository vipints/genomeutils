#!/usr/bin/env python 
"""
master script to execute the pipeline
"""

import os 
import sys 
import yaml 

from optparse import OptionParser

from fetch_remote_data import download_data as dld
from fetch_remote_data import prepare_data as ppd

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
    parser.add_option( "-a", "--genome_index", action="store_true", dest="genome_index", default=False, help="Create genome index to align the reads.")
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
        download_public_data(config_file)

    if options.genome_index:
        create_genome_index(config_file)



def create_genome_index(yaml_config):
    """
    wrapper for calling genome index function 
    """

    from signal_labels import experiment_details_db as expdb

    orgdb = expdb.experiment_db(yaml_config)

    for org_name, det in orgdb.items():
        
        ## TODO parallel submission with pygrid

        ppd.create_star_genome_index(det['fasta'], det['genome_index_dir'], det['gtf'], 1, det['read_length']-1)

        print 'Index generated for genome %s' % org_name
    


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
