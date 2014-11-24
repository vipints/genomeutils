#!/usr/bin/env python 
"""
generating standard genome indices. 

Usage: 
import prepare_data as pd 
.__doc__
.__doc__

Requirement:
    STAR aligner 
"""

import os 
import re
import sys 
import shutil
import subprocess 


def create_star_genome_index(fasta_file, out_dir, num_workers=1, genome_anno=None):
    """
    Creating STAR genome index with or without using genome annotation

    @args fasta_file: reference genome sequence file .fasta format 
    @type fasta_file: str 
    @args out_dir: genome index binary file storage place  
    @type out_dir: str 
    @args num_workers: number of threads to run (default value = 1)
    @type num_workers: int 
    @args genome_anno: genome annotation file (optional) 
    @type genome_anno: str 
    """

    cli = 'STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --runThreadN %d' % (out_dir, fasta_file, num_workers) 
    cli = 'STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --runThreadN %d --sjdbGTFfile %s' % (out_dir, fasta_file, num_workers, genome_anno) 
    cli = 'STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --runThreadN %d --sjdbGTFfile %s --sjdbGTFtagExonParentTranscript Parent' % (out_dir, fasta_file, num_workers, genome_anno) 
    process = subprocess.Popen(cli, shell=True) 
    process.wait()


if __name__=="__main__":
    print __doc__
    
"""
"""
