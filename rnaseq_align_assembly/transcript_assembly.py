#!/usr/bin/env python 
"""
Program to run transcriptome assembly using TransriptSkimmer and cufflinks program on read alignment data 

Requirement: 
    TransriptSkimmer - 
    cufflinks - 

"""

import os 
import subprocess


def run_trsk(org_name, gio_file, bam_file, res_path, gff_file="pred_genes.gff", reg_file="pred_region.bed"):
    """
    run TransriptSkimmer for provided genome

    """

    res_path = "%s/%s/trans_pred/" % (res_path, org_name)

    max_exon_length = 
    max_intron_length = 
    max_intergenic_region = 

    options="-maxel %d -ss -reglen 0.66 -maxic %d -minic 20 -maxin %d -mm 2 -exm 3 -indt 150 -exd 20 -tf 0.5 -inscf 3 -excut 3 -toff 100 -el 15" % (max_exon_length, max_intergenic_region, max_intron_length)

    os.chdir(res_path) 

    gio_file = 
    trsk_path = 

    cli = "%s./infer_genes -gio %s -bam %s -gff %s -reg %s %s" % (trsk_path, gio_file, bam_file, gff_file, reg_file, options)  
    sys.stdout.write('\trun %s \n' % cli)

    process = subprocess.Popen(cli, shell=True) 
    process.wait()

