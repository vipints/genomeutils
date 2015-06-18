#!/usr/bin/env python 
"""
Program to run different search for protein domain sequence to protein sequence 

Requirement: 
    blat - http://
    
"""

import os
import sys
import time
import tempfile
import subprocess


def BLATrun(infile, outfile, data_type, outform):
    """
    BLAT run

    @args infile: fasta file with different genome sequence
    @type infile: str 
    @args outfile: multiple sequence alignments are reported 
    @type outfile: str 
    @args data_type: DNA, PROTEIN 
    @type data_type: str 
    @args outform: CLUSTAL, PHYLIP, FASTA 
    @type outform: str 
    """

    try:
        subprocess.call(["blat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        exit("Please make sure that the `blat` binary is in your $PATH")


    cli = 'blat -INFILE=%s -OUTFILE=%s -OUTORDER=%s -TYPE=%s -OUTPUT=%s' % (infile, outfile, out_order, data_type, outform)

    try:
        outlog = "blat_run-%s.log" % time.strftime("%Y_%m_%d_%H-%M-%S")
        tlf = open(outlog,'w')

        process = subprocess.Popen(cli, shell=True, stderr=tlf, stdout=tlf)
        returncode = process.wait()

        if returncode !=0:
            raise Exception, "Exit status return code = %i" % returncode

        sys.stdout.write("blat run finished.\n")
        tlf.close() 

    except Exception, e:
        sys.stdout.write('Error running blat.\n%s\n' %  str( e ))
        sys.exit(0)

