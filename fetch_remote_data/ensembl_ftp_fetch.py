#!/usr/bin/env python 
"""
Download genome sequence from ENSEMBL ftp page.
"""
import sys 
import re, os 
import urllib2
import shutil

# update the ENSEMBL release version and downloading path  
ensembl_release_version = str(69)
download_path = "/tmp/"

base_url_fasta = "ftp://ftp.ensembl.org/pub/release-"+ ensembl_release_version +"/fasta/"

org_file = urllib2.urlopen(base_url_fasta)
for org_name in org_file:
    org_name=org_name.strip("\n\r")
    if org_name.split()[-1] in ['ailuropoda_melanoleuca', 'ancestral_alleles']:
        continue
    sub_path=org_name.split()[-1]+'/dna/'
    print org_name.split()[-1]
    fa_files = urllib2.urlopen(base_url_fasta+sub_path)
    for fa_name in fa_files:
        fa_name =fa_name.strip('\n\r')
        if re.search(r'.*.dna.toplevel.fa.gz$', fa_name.split()[-1]):
            os.makedirs(download_path+org_name.split()[-1]+"/ensembl_release-"+ensembl_release_version)
            tempfile=open(download_path+org_name.split()[-1]+"/ensembl_release-"+ensembl_release_version+"/"+fa_name.split()[-1], "wb")
            ftp_file=urllib2.urlopen(base_url_fasta+sub_path+fa_name.split()[-1])
            sys.stdout.write('\tdownloading '+fa_name.split()[-1]+' ... ')
            shutil.copyfileobj(ftp_file, tempfile)
            tempfile.close()
            ftp_file.close()
            sys.stdout.write("done\n")
            break
    fa_files.close()
org_file.close()
