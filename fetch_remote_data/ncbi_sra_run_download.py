#!/usr/bin/env python 
"""
Download RNA-sequncing reads trace file from NCBI SRA. 
This program will be able to download .SRA files based on the Run id. 

http://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR549/SRR549333/SRR549333.sra
or 
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR549/SRR549333/SRR549333.sra

"""
import urllib2
import shutil
import re, os

# TODO ERA url to be added, 
base_url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/" 


