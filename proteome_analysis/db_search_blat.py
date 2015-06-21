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


def BLATrun(database, query_seq, data_type, outform):
    """
    BLAT run

    @args database: fasta file with different sequence 
    @type database: str 
    @args query_seq: fasta file with query sequence to search 
    @type query_seq: str 
    @args data_type: DNA, PROTEIN, RNA 
    @type data_type: str 
    @args outform: PSL format filename 
    @type outform: str 
    """

    try:
        subprocess.call(["blat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        exit("Please make sure that the `blat` binary is in your $PATH")

    cli = 'blat %s %s -t=%s -q=%s %s' % (database, query_seq, data_type, data_type, outform)

    try:
        outlog = "blat_run-%s.log" % time.strftime("%Y_%m_%d_%H-%M-%S")
        tlf = open(outlog,'w')

        process = subprocess.Popen(cli, shell=True, stderr=tlf, stdout=tlf)
        returncode = process.wait()

        if returncode !=0:
            raise Exception, "Exit status return code = %i" % returncode

        sys.stdout.write("blat run finished. result saved at %s\n" % outform)
        tlf.close() 

    except Exception, e:
        sys.stdout.write('Error running blat.\n%s\n' %  str( e ))
        sys.exit(0)

