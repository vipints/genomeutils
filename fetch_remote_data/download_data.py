#!/usr/bin/env python 
"""
modules for downloading RNA-sequencing reads trace file from 
NCBI Short Read Archive, genome sequence and genome annotations 
from ENEMBL and Phytozome server. 

Usage: 
import download_data as dl 
dl.download_sra_file.__doc__
dl.uncompress_sra_file.__doc__

Requirement:
    fastq-dump - sratoolkit: http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software
"""

import os 
import re
import sys 
import shutil
import urllib2
import subprocess 

class MyException( Exception ): 
    pass

def fetch_phytozome_gff(release_version, species_name, download_path):
    """
    Download genome sequence from Phytozome ftp page.

    @args release_version: release version 
    @type release_version: str 
    @args species_name: organism name 
    @type species_name: str 
    @args download_path: file download path 
    @type download_path: str 
    """
    ## check the url for getting the recent version of the repository 
    base_url_gff = 'ftp://ftp.jgi-psf.org/pub/compgen/phytozome/%s/' % release_version
    
    try:
        org_names = urllib2.urlopen(base_url_gff)
    except urllib2.URLError, err_release:
        print "phytozome_release_version %s is NOT found" % release_version
        print err_release
        sys.exit(-1)

    for ORG in org_names:
        ORG = ORG.strip("\n\r")

        if ORG.split()[-1] != species_name:
            continue

        ## updating the base_url 
        base_url_gff = '%s%s/annotation/' % (base_url_gff, species_name)
        try:
            gff_files = urllib2.urlopen(base_url_gff)
        except urllib2.URLError, err_faseq:
            print "phytozome_release genome annotation missing" % base_url_gff
            print err_faseq
            sys.exit(-1)

        ## mapping to short names  Athaliana --> A_thaliana
        org_short_name = "%s_%s" % (species_name[0], species_name[1:])

        ## setting up the download path 
        ## download the files in ex: /home/tmp/V_vinifera/phytozome_v9.0/Vvinifera_145.fa.gz
        base_file_path = "%s/%s/phytozome_%s" % (download_path, org_short_name, release_version)
        if not os.path.exists(base_file_path):
            os.makedirs(base_file_path)

        for gff_name in gff_files:
            gff_name =gff_name.strip('\n\r')

            if re.search(r'.*_\d+_gene.gff3.gz$', gff_name.split()[-1]):
                tempfile=open("%s/%s" % (base_file_path, gff_name.split()[-1]), "wb")
                
                try:
                    ftp_file=urllib2.urlopen(base_url_gff+gff_name.split()[-1])
                except urllib2.URLError, err_file:
                    print err_file
                    sys.exit(-1)

                sys.stdout.write('\tdownloading %s ... ' % gff_name.split()[-1])
                shutil.copyfileobj(ftp_file, tempfile)

                tempfile.close()
                ftp_file.close()

                sys.stdout.write("done\n")
        gff_files.close()
    org_names.close()


def fetch_phytozome_fasta(release_version, species_name, download_path):
    """
    Download genome sequence from Phytozome ftp page.

    @args release_version: release version 
    @type release_version: str 
    @args species_name: organism name 
    @type species_name: str 
    @args download_path: file download path 
    @type download_path: str 
    """
    ## check the url for getting the recent version of the repository 
    base_url_fasta = 'ftp://ftp.jgi-psf.org/pub/compgen/phytozome/%s/' % release_version
    
    try:
        org_names = urllib2.urlopen(base_url_fasta)
    except urllib2.URLError, err_release:
        print "phytozome_release_version %s is NOT found" % release_version
        print err_release
        sys.exit(-1)

    for ORG in org_names:
        ORG = ORG.strip("\n\r")

        if ORG.split()[-1] != species_name:
            continue

        ## updating the base_url 
        base_url_fasta = '%s%s/assembly/' % (base_url_fasta, species_name)
        try:
            fa_files = urllib2.urlopen(base_url_fasta)
        except urllib2.URLError, err_faseq:
            print "phytozome_release genome sequence missing" % base_url_fasta
            print err_faseq
            sys.exit(-1)

        ## mapping to short names  Athaliana --> A_thaliana
        org_short_name = "%s_%s" % (species_name[0], species_name[1:])

        ## setting up the download path 
        ## download the files in ex: /home/tmp/V_vinifera/phytozome_v9.0/Vvinifera_145.fa.gz
        base_file_path = "%s/%s/phytozome_%s" % (download_path, org_short_name, release_version)
        if not os.path.exists(base_file_path):
            os.makedirs(base_file_path)

        for fa_name in fa_files:
            fa_name =fa_name.strip('\n\r')

            if re.search(r'.*_\d+.fa.gz$', fa_name.split()[-1]):
                tempfile=open("%s/%s" % (base_file_path, fa_name.split()[-1]), "wb")
                
                try:
                    ftp_file=urllib2.urlopen(base_url_fasta+fa_name.split()[-1])
                except urllib2.URLError, err_file:
                    print err_file
                    sys.exit(-1)

                sys.stdout.write('\tdownloading %s ... ' % fa_name.split()[-1])
                shutil.copyfileobj(ftp_file, tempfile)

                tempfile.close()
                ftp_file.close()

                sys.stdout.write("done\n")
        fa_files.close()
    org_names.close()

def fetch_ensembl_metazoa_gtf(release_version, species_name, download_path):
    """
    Download genome annotation from ENSEMBLGENOMES ftp page.
    
    @args release_version: ensembl release version 
    @type release_version: str 
    @args species_name: organism name (example: anopheles_gambiae)
    @type species_name: str 
    @args download_path: file download path 
    @type download_path: str 
    """

    ## check the url for getting the recent version of the repository 
    base_url_gtf = "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-%s/gtf/" % release_version 
    
    try:
        org_file = urllib2.urlopen(base_url_gtf)
    except urllib2.URLError, err_release:
        print "ensembl_release_version %s is NOT found" % release_version
        print err_release
        sys.exit(-1)

    for org_name in org_file:
        org_name=org_name.strip("\n\r")

        ## check the organism directory at ftp remote folder 
        if org_name.split()[-1] != species_name: 
            continue

        ## updating the base url 
        base_url_gtf = '%s%s/' % (base_url_gtf, org_name.split()[-1])

        try:
            gtf_files = urllib2.urlopen(base_url_gtf)
        except urllib2.URLError, err_gtf:
            print "ensembl_release genome annotation missing" % base_url_gtf
            print err_gtf
            sys.exit(-1)

        ## mapping to short names  Arabidopsis_thaliana --> A_thaliana
        genus, species = species_name.strip().split("_")
        org_short_name = "%s_%s" % (genus[0].upper(), species)

        base_file_path = "%s/%s/ensembl_release_%s" % (download_path, org_short_name, release_version)
        if not os.path.exists(base_file_path):
            os.makedirs(base_file_path)
            
        for gtf_name in gtf_files:
            gtf_name =gtf_name.strip('\n\r')

            if re.search(r'.*.\d+.gtf.gz$', gtf_name.split()[-1]):
                tempfile=open("%s/%s" % (base_file_path, gtf_name.split()[-1]), "wb")

                try:
                    ftp_file=urllib2.urlopen(base_url_gtf+gtf_name.split()[-1])
                except urllib2.URLError, err_file:
                    print err_file
                    sys.exit(-1)

                sys.stdout.write('\tdownloading %s ...\n' % gtf_name.split()[-1])
                shutil.copyfileobj(ftp_file, tempfile)
                sys.stdout.write('\tsaved at %s/%s \n' % (base_file_path, gtf_name.split()[-1]))

                tempfile.close()
                ftp_file.close()

        gtf_files.close()
    org_file.close()


def fetch_ensembl_gtf(release_version, species_name, download_path):
    """
    Download genome annotation from ENSEMBL ftp page.
    
    @args release_version: ensembl release version 
    @type release_version: str 
    @args species_name: organism name (example: homo_sapiens)
    @type species_name: str 
    @args download_path: file download path 
    @type download_path: str 
    """
    ## check the url for getting the recent version of the repository 
    base_url_gtf = "ftp://ftp.ensembl.org/pub/release-%s/gtf/" % release_version 

    try:
        org_file = urllib2.urlopen(base_url_gtf)
    except urllib2.URLError, err_release:
        print "ensembl_release_version %s is NOT found" % release_version
        print err_release
        sys.exit(-1)

    for org_name in org_file:
        org_name=org_name.strip("\n\r")
        
        ## check the organism directory at ftp remote folder 
        if org_name.split()[-1] != species_name: 
            continue

        ## updating the base url 
        base_url_gtf = '%s%s/' % (base_url_gtf, org_name.split()[-1])
        
        try:
            gtf_files = urllib2.urlopen(base_url_gtf)
        except urllib2.URLError, err_gtf:
            print "ensembl_release genome annotation missing" % base_url_gtf
            print err_gtf
            sys.exit(-1)

        ## mapping to short names  Arabidopsis_thaliana --> A_thaliana
        genus, species = species_name.strip().split("_")
        org_short_name = "%s_%s" % (genus[0].upper(), species)

        ## setting up the download path 
        base_file_path = "%s/%s/ensembl_release_%s" % (download_path, org_short_name, release_version)
        if not os.path.exists(base_file_path):
            os.makedirs(base_file_path)

        for gtf_name in gtf_files:
            gtf_name =gtf_name.strip('\n\r')

            if re.search(r'.*.\d+.gtf.gz$', gtf_name.split()[-1]):
                tempfile=open("%s/%s" % (base_file_path, gtf_name.split()[-1]), "wb")

                try:
                    ftp_file=urllib2.urlopen(base_url_gtf+gtf_name.split()[-1])
                except urllib2.URLError, err_file:
                    print err_file
                    sys.exit(-1)

                sys.stdout.write('\tdownloading %s ...\n' % gtf_name.split()[-1])
                shutil.copyfileobj(ftp_file, tempfile)
                sys.stdout.write('\tsaved at %s/%s \n' % (base_file_path, gtf_name.split()[-1]))

                tempfile.close()
                ftp_file.close()

        gtf_files.close()
    org_file.close()


def fetch_ensembl_metazoa_fasta(release_version, species_name, download_path):
    """
    Download genome sequence from ensemblgenomes ftp page.
    
    @args release_version: ensembl release version 
    @type release_version: str 
    @args species_name: organism name (example: homo_sapiens)
    @type species_name: str 
    @args download_path: file download path 
    @type download_path: str 
    """

    ## check the url for getting the recent version of the repository 
    base_url_fasta = "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-%s/fasta/" % release_version 

    try:
        org_file = urllib2.urlopen(base_url_fasta)
    except urllib2.URLError, err_release:
        print "ensembl_release_version %s is NOT found" % release_version
        print err_release
        sys.exit(-1)

    for org_name in org_file:
        org_name=org_name.strip("\n\r")
        
        ## check the organism directory at ftp remote folder 
        if org_name.split()[-1] != species_name: 
            continue

        ## updating the base url 
        base_url_fasta = '%s%s/dna/' % (base_url_fasta, org_name.split()[-1])

        try:
            fa_files = urllib2.urlopen(base_url_fasta)
        except urllib2.URLError, err_faseq:
            print "ensembl_release genome sequence missing" % base_url_fasta
            print err_faseq
            sys.exit(-1)

        ## mapping to short names  Arabidopsis_thaliana --> A_thaliana
        genus, species = species_name.strip().split("_")
        org_short_name = "%s_%s" % (genus[0].upper(), species)

        base_file_path = "%s/%s/ensembl_release_%s" % (download_path, org_short_name, release_version)
        if not os.path.exists(base_file_path):
            os.makedirs(base_file_path)

        for fa_name in fa_files:
            fa_name =fa_name.strip('\n\r')

            ## include repeatmasked genome 
            if re.search(r'.*.dna.toplevel.fa.gz$', fa_name.split()[-1]) or re.search(r'.*.dna_rm.toplevel.fa.gz$', fa_name.split()[-1]):
                tempfile=open("%s/%s" % (base_file_path, fa_name.split()[-1]), "wb")

                try:
                    ftp_file=urllib2.urlopen(base_url_fasta+fa_name.split()[-1])
                except urllib2.URLError, err_file:
                    print err_file
                    sys.exit(-1)

                sys.stdout.write('\tdownloading %s ... \n' % fa_name.split()[-1])
                shutil.copyfileobj(ftp_file, tempfile)
                sys.stdout.write('\tsaved at %s/%s \n' % (base_file_path, fa_name.split()[-1]))

                tempfile.close()
                ftp_file.close()

        fa_files.close()
    org_file.close()


def fetch_ensembl_fasta(ensembl_release_version, species_name, download_path):
    """
    Download genome sequence from ENSEMBL ftp page.
    
    @args ensembl_release_version: ensembl release version 
    @type ensembl_release_version: str 
    @args species_name: organism name (example: homo_sapiens)
    @type species_name: str 
    @args download_path: file download path 
    @type download_path: str 
    """
    ## check the url for getting the recent version of the repository 
    base_url_fasta = "ftp://ftp.ensembl.org/pub/release-%s/fasta/" % ensembl_release_version 

    try:
        org_file = urllib2.urlopen(base_url_fasta)
    except urllib2.URLError, err_release:
        print "ensembl_release_version %s is NOT found" % ensembl_release_version
        print err_release
        sys.exit(-1)

    for org_name in org_file:
        org_name=org_name.strip("\n\r")
        
        ## check the organism directory at ftp remote folder 
        if org_name.split()[-1] != species_name: 
            continue

        ## updating the base url 
        base_url_fasta = '%s%s/dna/' % (base_url_fasta, org_name.split()[-1])
        
        try:
            fa_files = urllib2.urlopen(base_url_fasta)
        except urllib2.URLError, err_faseq:
            print "ensembl_release genome sequence missing" % base_url_fasta
            print err_faseq
            sys.exit(-1)

        ## mapping to short names  Arabidopsis_thaliana --> A_thaliana
        genus, species = species_name.strip().split("_")
        org_short_name = "%s_%s" % (genus[0].upper(), species)

        ## download the files in ex: /home/tmp/F_albicollis/ensembl_release-77/Ficedula_albicollis.FicAlb_1.4.dna_rm.toplevel.fa.gz
        base_file_path = "%s/%s/ensembl_release_%s" % (download_path, org_short_name, ensembl_release_version)
        if not os.path.exists(base_file_path):
            os.makedirs(base_file_path)

        for fa_name in fa_files:
            fa_name =fa_name.strip('\n\r')

            ## include repeatmasked genome 
            if re.search(r'.*.dna.toplevel.fa.gz$', fa_name.split()[-1]) or re.search(r'.*.dna_rm.toplevel.fa.gz$', fa_name.split()[-1]):
                tempfile=open("%s/%s" % (base_file_path, fa_name.split()[-1]), "wb")

                try:
                    ftp_file=urllib2.urlopen(base_url_fasta+fa_name.split()[-1])
                except urllib2.URLError, err_file:
                    print err_file
                    sys.exit(-1)

                sys.stdout.write('\tdownloading %s ... ' % fa_name.split()[-1])
                shutil.copyfileobj(ftp_file, tempfile)

                tempfile.close()
                ftp_file.close()

                sys.stdout.write("done\n")
        fa_files.close()
    org_file.close()


def download_sra_file(RUNID, download_path):
    """
    Download the SRA file

    @args RUNID: SRA run ID 
    @type RUNID: str 
    @args download_path: SRA file download path 
    @type download_path: str 
    """
    ## ncbi sra trace url 
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

    return out_file_name


def uncompress_sra_file(out_file_name, download_path, lib_type="pe", out_compress="bzip2"):
    """
    Uncompress downloaded SRA file 

    @args out_file_name: Downloaded SRA file name 
    @type out_file_name: str 
    @args download_path: SRA file uncompressing path 
    @type download_path: str 
    @args lib_type: Library layout (default pe) 
    @type lib_type: str 
    @args out_compress: compress format for result file (default bzip2)
    @type out_compress: str 

    NOTE: This module expects sratoolkit is available under PATH variable or add the path below line.
    """
    ## add the installation path of sratoolkit  
    #os.environ['PATH'] += os.pathsep + '/home/share/software/sratoolkit/sratoolkit.2.3.1-centos_linux64/bin/'

    ## depends on the compress type and library protocol type
    if lib_type in ['pe', 'PE', 'paired-end']:
        cli = 'fastq-dump \
        --%s \
        --split-3 \
        --outdir %s \
        %s' % (out_compress, download_path, out_file_name) 
    elif lib_type in ['se', 'SE', 'single-end']:
        cli = 'fastq-dump \
        --%s \
        --outdir %s \
        %s' % (out_compress, download_path, out_file_name)
    else:
        print 'Error! Library layout [PE/pe/paired-end|SE/se/single-end]'
        print '\tYour Library layout is ', lib_type
        sys.exit(-1)

    ## split the .SRA format file based on the library layout
    sys.stdout.write('\trun %s \n' % cli)
    process = subprocess.Popen(cli, shell=True) 
    process.wait()


if __name__=="__main__":
    print __doc__
    
"""
    sra_run_id = ""
    sra_data_path = ""
    library_type = ""
    compress_format = ""

    sra_file = download_sra_file(sra_run_id, sra_data_path)
    uncompress_sra_file(sra_file, data_path, library_type, compress_format)

    release_num = ""
    org_name = "" # homo_sapiens, Vcarteri
    fasta_data_path = ""
    
    fetch_ensembl_fasta(release_num, org_name, fasta_data_path)
    fetch_phytozome_fasta(release_num, org_name, fasta_data_path)

    fetch_phytozome_gff(release_num, org_name, fasta_data_path)
    fetch_ensembl_gtf(release_num, org_name, fasta_data_path)
"""
