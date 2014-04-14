#!/usr/bin/env python 
"""
Program to run clustalw2 

Usage: python clustalw.py infile outfile data_type out_format 

INFILE: FASTA
DATATYPE: DNA, PROTEIN 
OUTFORM: CLUSTAL, PHYLIP, FASTA

Requirement: 
    install clustalw2 and set the path correctly 

TODO: include more commandline features to clustalw2
check the inputs from user 
"""

import os
import sys
import tempfile
import subprocess

def ClustalWrun(infile, outfile, data_type, outform):
    """
    clustalw2 run 
    """
    
    outlog = "clustalw2_run.log"
    tlf = open(outlog,'w')

    out_order = "INPUT" # ALIGNED

    cl = ['clustalw2 -INFILE=%s -OUTFILE=%s -OUTORDER=%s -TYPE=%s -OUTPUT=%s' % (infile, outfile, out_order, data_type, outform)]
    process = subprocess.Popen(' '.join(cl), shell=True, stderr=tlf, stdout=tlf)
    rval = process.wait()
    tlf.close() 


try:
    INFILE = sys.argv[1]
    OUTFILE = sys.argv[2]
    DATATYPE = sys.argv[3]
    OUTFORM = sys.argv[4]
except:
    print __doc__ 
    sys.exit(-1) 



ClustalWrun(INFILE, OUTFILE, DATATYPE, OUTFORM) 

