#!/usr/bin/env python 
"""
Download RNA-sequencing reads trace file from NCBI Short Read Archive repository based on Run ID. 

Usage: python ncbi_sra_run_download.py run_id out_dir pe/se{paired-end/single-end}  

Requirement:
    fastq-dump - sratoolkit: http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software
"""

import os 
import sys 
import shutil
import urllib2
import subprocess 

class MyException( Exception ): 
    pass


def download_sra_file(RUNID=None, download_path=None, lib_type='pe', out_compress="bzip2"):
    """
    Download the SRA file

    @args RUNID: SRA run ID 
    @type RUNID: str 
    @args download_path: SRA file download path 
    @type download_path: str 
    @args lib_type: Library layout 
    @type lib_type: str 
    @args out_compress: compress format for result file  
    @type out_compress: str 
    """

    base_url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/"

    ## sanity check for run id 
    assert len(RUNID) in [9, 10], 'Error in SRA Run ID format [ex: SRR548309, SRR1050788] -- %s --' % RUNID

    ## build the complete url based on the RunID 
    if not RUNID[0:3] in ['DRR', 'SRR', 'ERR']:
        print 'Error! Experiment Run ID must start with DRR or SRR or ERR'
        print '\tYour Run ID start with ', RUNID[0:3], ' prefix'
        print '\tProgram cannot continue, Exiting...'
        sys.exit(-1)

    ## adding sub folder 
    base_url = '%s%s/' % (base_url, RUNID[0:3]) 

    try:
        isinstance(int(RUNID[3:6]), int)
    except:
        print 'Error! Experiment Run ID will be in the format {SRR/ERR/DRR}NNN'
        print '\tYour Run ID is ', RUNID[0:6]
        print '\tProgram cannot continue, Exiting...'
        sys.exit(-1)
    ## adding sub - sub folder - ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR105/ 
    base_url = '%s%s/' % (base_url, RUNID[0:6]) 

    try:
        isinstance(int(RUNID[6:10]), int)
        ## adding sub - sub sub folder ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR105/SRR1050788/SRR1050788.sra
        base_url = '%s%s/%s.sra' % (base_url, RUNID[0:10], RUNID[0:10])
   
    except:
        try:
            isinstance(int(RUNID[6:9]), int)
        except:
            print 'Error! Experiment Run ID will be in the format {SRR/ERR/DRR}NNNNNN'
            print '\tYour Run ID is ', RUNID[0:9]
            print '\tProgram cannot continue, Exiting...'
            sys.exit(-1)

    ## adding sub - sub sub folder - ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR548/SRR548309/SRR548309.sra 
    base_url = '%s%s/%s.sra' % (base_url, RUNID[0:9], RUNID[0:9]) 

    ## create downloadpath if doesnot exists 
    if not os.path.exists(download_path):
        os.makedirs(download_path)

    ## getting the absolute path 
    download_path = os.path.realpath(download_path)

    ## out file name 
    out_file_name = '%s/%s.sra' % (download_path, RUNID)

    ## creating temporory file 
    tempfile=open(out_file_name, "wb")

    ## start fetching the remote file 
    try:
        sra_file = urllib2.urlopen(base_url, timeout=1 )
    except urllib2.URLError, err:
        raise MyException("There is an Error: %r" % err)

    ## download job starts 
    sys.stdout.write('\tdownloading %s ...\n' % base_url)
    shutil.copyfileobj(sra_file, tempfile)
    sys.stdout.write('\tsaved at %s done!\n' % out_file_name)

    tempfile.close()
    sra_file.close()


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
sys.stdout.write('\trun %s \n' % cli)
process = subprocess.Popen(cli, shell=True) 
process.wait()


if __name__=="__main__":
    
    try:
        RUNID = sys.argv[1]
        download_path = sys.argv[2]
        lib_type = sys.argv[3]
    except:
        print __doc__
        sys.exit(-1)
