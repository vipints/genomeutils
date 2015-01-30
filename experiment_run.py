#!/usr/bin/env python 
"""
master script to execute the pipeline
"""

import os 
import sys 

from signal_labels import org_details_db as odb
from fetch_remote_data import download_data as dld


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



if __name__=="__main":
    
    infile = ""
    data_path = ""
    exp_path = ""

    org_details = odb.make_org_db(infile, data_path, exp_path) 

    download_fasta(org_details) 

    download_gtf(org_details)     
   

