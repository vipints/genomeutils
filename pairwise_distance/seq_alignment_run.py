#!/usr/bin/env python 
"""
Program to run different multiple sequence alignment programs. 

Requirement: 
    mafft - http://mafft.cbrc.jp/alignment/server/index.html
    clustalw2 - http://www.ebi.ac.uk/Tools/msa/clustalw2/ 

    install these packages and set the path correctly if it is not a standard one. 
"""

import os
import sys
import time
import tempfile
import subprocess


def MAFFTrun(infile, outfile, threads=2):
    """
    mafft run command line 
    
    @args infile: fasta file with different genome sequence
    @type infile: str 
    @args outfile: multiple sequence alignments are reported, example: CLUSTAL, PHYLIP, FASTA
    @type outfile: str 
    @args threads: number of cores used for executing the program
    @type threads: int 
    """
    
    outlog = "mafft_run-%s.log" % time.strftime("%Y_%m_%d_%H-%M-%S") 
    tlf = open(outlog, 'w')

    # TODO: include more commandline features to mafft 
    out_order = "ALIGNED" # "INPUT" # ALIGNED

    cl = ['mafft --thread %d --threadit %d --reorder --anysymbol --auto -OUTPUT=%s %s' % (threads, 0, outfile, infile)]
    process = subprocess.Popen(' '.join(cl), shell=True, stderr=tlf, stdout=tlf)
    rval = process.wait()
    tlf.close() 


def ClustalWrun(infile, outfile, data_type, outform):
    """
    clustalw2 run 

    @args infile: fasta file with different genome sequence
    @type infile: str 
    @args outfile: multiple sequence alignments are reported 
    @type outfile: str 
    @args data_type: DNA, PROTEIN 
    @type data_type: str 
    @args outform: CLUSTAL, PHYLIP, FASTA 
    @type outform: str 
    """
    
    outlog = "clustalw2_run-%s.log" % time.strftime("%Y_%m_%d_%H-%M-%S")
    tlf = open(outlog,'w')
    
    #TODO: include more commandline features to clustalw2
    out_order = "ALIGNED" # "INPUT" # ALIGNED

    cl = ['clustalw2 -INFILE=%s -OUTFILE=%s -OUTORDER=%s -TYPE=%s -OUTPUT=%s' % (infile, outfile, out_order, data_type, outform)]
    process = subprocess.Popen(' '.join(cl), shell=True, stderr=tlf, stdout=tlf)
    rval = process.wait()
    tlf.close() 

