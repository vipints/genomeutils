#!/usr/bin/env python 
"""
Program to run MAFFT 

Usage: python mafft_run.py infile outfile data_type out_format 

INFILE: FASTA

Requirement: 
    install mafft and set the path correctly 

TODO: include more commandline features to mafft 
check the inputs from user 
http://mafft.cbrc.jp/alignment/server/index.html
"""

import os
import sys
import tempfile
import subprocess

def MAFFTrun(infile, outfile):
    """
    mafft run 
    """
    
    outlog = "mafft_run.log"
    tlf = open(outlog,'w')

    out_order = "ALIGNED" # "INPUT" # ALIGNED

    cl = ['mafft --thread %d --threadit %d --reorder --anysymbol --auto -OUTPUT=%s %s' % (2, 0, outfile, infile)]
    process = subprocess.Popen(' '.join(cl), shell=True, stderr=tlf, stdout=tlf)
    rval = process.wait()
    tlf.close() 


try:
    INFILE = sys.argv[1]
    OUTFILE = sys.argv[2]
except:
    print __doc__ 
    sys.exit(-1) 


MAFFTrun(INFILE, OUTFILE) 

