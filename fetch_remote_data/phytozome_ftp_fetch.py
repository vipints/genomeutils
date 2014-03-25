#!/usr/bin/env python 
"""
Download genome sequence from Phytozome ftp page.
"""
import sys 
import re, os 
import urllib2
import shutil

# update the release version and downloading path  
release_version = 'v9.0'
download_path = "/tmp/"

base_url = 'ftp://ftp.jgi-psf.org/pub/compgen/phytozome/%s' % release_version

org_names = urllib2.urlopen(base_url)
for org_name in org_names:
    org_name = org_name.strip("\n\r")
    print org_name.split()[-1]

org_file.close()
