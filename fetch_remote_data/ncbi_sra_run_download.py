#!/usr/bin/env python 
"""
Download RNA-sequencing reads trace file from NCBI Short Read Archive repository based on Run ID. 

Usage: python ncbi_sra_run_download.py run_id out_dir pe/se{paired-end/single-end}  

Requirement:
    fastq-dump - sratoolkit: http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software
"""


import re 
import os 
import sys 
import shutil
import urllib2
import subprocess 

class MyException( Exception ): 
    pass

base_url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/"

try:
    RUNID = sys.argv[1]
    download_path = sys.argv[2]
    lib_type = sys.argv[3]
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

## getting the absolute path 
download_path = os.path.realpath(download_path)

## out file name 
out_file_name = download_path + '/' + RUNID[0:9] + '.sra' 

## creating temporory file 
tempfile=open(out_file_name, "wb")

## start fetching the remote file 
try:
    sra_file = urllib2.urlopen(base_url, timeout=1 )
except urllib2.URLError, err:
    raise MyException("There is an Error: %r" % err)

## download job starts 
sys.stdout.write('\tdownloading ' + base_url + ' ... \n')
shutil.copyfileobj(sra_file, tempfile)
sys.stdout.write('\tsaved at ' + out_file_name + ' done! \n')

tempfile.close()
sra_file.close()

## default compression type
out_compress = "bzip2"

## depends on the compress type and library protocol type
if lib_type in ['pe', 'PE', 'paired-end']:
    cli = 'fastq-dump --%s --split-3 --outdir %s %s' % (out_compress, download_path, out_file_name) 
elif lib_type in ['se', 'SE', 'single-end']:
    cli = 'fastq-dump --%s --outdir %s %s' % (out_compress, download_path, out_file_name)
else:
    print 'Error! Library layout [PE/pe/paired-end|SE/se/single-end]'
    print '\tYour Library layout is ', lib_type
    print '\tProgram cannot continue, Exiting...'
    sys.exit(-1)

## add the installation path of sratoolkit  
os.environ['PATH'] += os.pathsep + '/home/share/software/sratoolkit/sratoolkit.2.3.1-centos_linux64/bin/'

## split the .SRA format file based on the library layout
sys.stdout.write('\trun ' + cli + '\n')
process = subprocess.Popen(cli, shell=True) 
process.wait()
