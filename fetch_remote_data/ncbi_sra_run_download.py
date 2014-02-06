#!/usr/bin/env python 
"""
Download RNA-sequencing reads trace file from NCBI Short Read Archive repository based on Run ID. 

Usage: python ncbi_sra_run_download.py run_id out_dir pe/se{paired-end/single-end}  
"""

import sys 
import re, os
import urllib2
import shutil

base_url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/"

try:
    RUNID = sys.argv[1]
    download_path = sys.argv[2]
except:
    print __doc__
    sys.exit(-1)

## sanity check for run id 
assert len(RUNID)==9, 'Error in SRA Run ID format [ex: SRR548309] '+ RUNID

## build the complete url based on the RunID 
if not RUNID[0:3] in ['DRR', 'SRR', 'ERR']:
    print 'Error! Experiment Run ID must start with DRR or SRR or ERR'
    print '\tYour Run ID start with ', RUNID[0:3], ' prefix'
    print '\tProgram cannot continue, Exiting...'
    sys.exit(-1)

## adding sub folder 
base_url = base_url + RUNID[0:3] + '/'

try:
    isinstance(int(RUNID[3:6]), int)
except:
    print 'Error! Experiment Run ID will be in the format {SRR/ERR/DRR}NNN'
    print '\tYour Run ID is ', RUNID[0:6]
    print '\tProgram cannot continue, Exiting...'
    sys.exit(-1)

## adding sub - sub folder 
base_url = base_url + RUNID[0:6] + '/'

try:
    isinstance(int(RUNID[6:9]), int)
except:
    print 'Error! Experiment Run ID will be in the format {SRR/ERR/DRR}NNNNNN'
    print '\tYour Run ID is ', RUNID[0:9]
    print '\tProgram cannot continue, Exiting...'
    sys.exit(-1)

## adding sub - sub sub folder 
base_url = base_url + RUNID[0:9] + '/'

## append the .sra extension to the run id 
base_url = base_url + RUNID[0:9] + '.sra'

## create downloadpath if doesnot exists 
if not os.path.exists(download_path):
    os.makedirs(download_path)

## creating temporory file 
tempfile=open(download_path + "/" + RUNID[0:9] + '.sra', "wb")

## start fetching the remote file 
sra_file = urllib2.urlopen(base_url)

## download job starts 
sys.stdout.write('\tdownloading ' + base_url + ' ... \n')
shutil.copyfileobj(sra_file, tempfile)
sys.stdout.write('\tsaved at ' + download_path + '/' + RUNID[0:9] + '.sra done! \n')

tempfile.close()
sra_file.close()


