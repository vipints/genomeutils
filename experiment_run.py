#!/usr/bin/env python 
"""
master script to control different operations in training examples generating pipeline.

usage:
    python experiment_run.py <YAML config> -h 

    sample yaml config file located at config/

requirement:
    pygridtools for distributed computing 
    packages/modules depends on the operation  
"""

import os 
import sys 
import yaml 

try:
    import libpyjobrunner as pg
except:
    sys.stdout.write('warning: pygridtools are not available, distributed computing task will be disrupted\n')

from optparse import OptionParser
from signal_labels import experiment_details_db as expdb

assert sys.version_info[:2] >= ( 2, 4 )

def main():
    """
    Managing the experiment run in different operation mode: 

    Options

    -1 download_sra file from NCBI SRA service
    -d decompose_sra decompress the SRA file  
    -a annotation fetch genome annotation file from public database servers mainly ensembl, phytozome 
    -g genome fetch genome sequence file from public database

    manual cleaning of genome sequence and annotation 
    
    -2 genome_index create genome indices for STAR alignment 
    -i insert_size calculate the insert size based on the raw fastq files  
    -3 read_mapping aligning reads to the genome using STAR aligner 
    -m multi_map resolving the right place for read mapped in multiple location on the genome
    -u uniq_read recover uniquely aligned reads from the star alignment 
    -4 trsk_pred transcript prediction using TranscriptSkimmer
    -c cuff_pred transcript assembly by Cufflinks 
    -s stringtie_pred transcript assembly by StringTie
    -5 filter_trsk applying filter to the trsk predicted gene models 
    -b filter_cuff applying filter to the cufflinks predicted gene models 
    -f filter_db applying filter to the online db genome annotations
    -6 trsk_label generating labels for genomic signal based on the trsk feature annotation  
    -t cuff_label generating labels for genomic signal based on the cufflinks feature annotation 
    -p db_label generating labels for genomic signal based on online db annotations
    """

    parser = OptionParser(usage='usage: %prog <YAML config> [required option]') 

    parser.add_option( "-1", "--download_sra", action="store_true", dest="download_sra", default=False, help="Download sra file based on run id from NCBI SRA/ENA repositories.")
    parser.add_option( "-d", "--decompose_sra", action="store_true", dest="decompose_sra", default=False, help="Decompress the sra file according to the library type.")
    parser.add_option( "-a", "--annotation", action="store_true", dest="annotation", default=False, help="Download genome annotation from public database resources.")
    parser.add_option( "-g", "--genome", action="store_true", dest="genome", default=False, help="Download genome sequence from public database resources.")
    parser.add_option( "-2", "--genome_index", action="store_true", dest="genome_index", default=False, help="Create STAR genome index based on genome sequence and annotations." )
    parser.add_option( "-i", "--insert_size", action="store_true", dest="insert_size", default=False, help="Calculate the library insert size from fastq files.")
    parser.add_option( "-3", "--read_mapping", action="store_true", dest="read_mapping", default=False, help="RNASeq read mapping to genome using STAR aligner." )
    parser.add_option( "-m", "--multi_map", action="store_true", dest="multi_map", default=False, help="MMR on aligned reads to resolve multimapping of reads." )
    parser.add_option( "-u", "--uniq_read", action="store_true", dest="uniq_read", default=False, help="Fetching uniquely mapped reads from bam file." )
    parser.add_option( "-4", "--trsk_pred", action="store_true", dest="trsk_pred", default=False, help="Transcript prediction using TranscriptSkimmer." )
    parser.add_option( "-c", "--cuff_pred", action="store_true", dest="cuff_pred", default=False, help="Transcript assembly using Cufflinks." )
    parser.add_option( "-s", "--stringtie_pred", action="store_true", dest="stringtie_pred", default=False, help="Transcript assembly using StringTie." )
    parser.add_option( "-5", "--filter_trsk", action="store_true", dest="filter_trsk", default=False, help="Apply filters to trsk predicted gene models." )
    parser.add_option( "-b", "--filter_cuff", action="store_true", dest="filter_cuff", default=False, help="Apply filter to the cufflinks predicted gene models." )
    parser.add_option( "-f", "--filter_db", action="store_true", dest="filter_db", default=False, help="Apply filter to the online db annotation gene models." )
    parser.add_option( "-6", "--trsk_label", action="store_true", dest="trsk_label", default=False, help="Fetch label sequences from TranscriptSkimmer annotations." )
    parser.add_option( "-t", "--cuff_label", action="store_true", dest="cuff_label", default=False, help="Fetch label sequences from cufflinks annotations." )
    parser.add_option( "-p", "--db_label", action="store_true", dest="db_label", default=False, help="Fetch label sequences from public online db annotation files." )

    ( options, args ) = parser.parse_args()
    try:
        config_file = args[0]
    except:
        exit(__doc__)

    if not (options.download_sra ^ options.decompose_sra ^ options.annotation ^ \
            options.genome ^ options.genome_index ^ options.insert_size ^ \
            options.read_mapping ^ options.multi_map ^ options.uniq_read ^ \
            options.trsk_pred ^ options.cuff_pred ^ options.filter_trsk ^ \
            options.trsk_label ^ options.filter_cuff ^ options.filter_db ^ \
            options.cuff_label ^ options.db_label ^ options.stringtie_pred):
        parser.print_help()
        sys.exit(-1)
        
    print('Using config file %s for the experiment.' % config_file)

    if options.download_sra:
        print 'Operation selected: Download sequencing reads file from ncbi-sra'
        download_sra_data(config_file)
    elif options.decompose_sra:
        print 'Operation selected: Decompress sra file'
        decompose_sra_file(config_file)
    elif options.annotation:
        print 'Operation selected: Downloading genome annotation file'
        download_gtf(config_file)
    elif options.genome:
        print 'Operation selected: Downloading genome sequence file'
        download_fasta(config_file)
    elif options.genome_index:
        print 'Operation selected: Create STAR genome index'
        create_genome_index(config_file)
    elif options.insert_size:
        print 'Operation selected: Calculate the library insert size from sequencing \
            read files'
        calculate_insert_size(config_file)
    elif options.read_mapping: 
        print 'Operation selected: Read alignment with STAR'
        align_rnaseq_reads(config_file)
    elif options.multi_map:
        print 'Operation selected: Multiple read mapper resolution with MMR'
        alignment_filter(config_file) 
    elif options.uniq_read:
        print 'Operation selected: Find uniquely mapped reads from star alignment'
        find_uniq_reads(config_file) 
    elif options.trsk_pred:
        print 'Operation selected: Transcript prediction based on mapped RNASeq read \
            data with TranscriptSkimmer'
        transcript_prediction_trsk(config_file)
    elif options.cuff_pred:
        print 'Operation selected: Transcript assembly based on mapped RNASeq read data \
            with Cufflinks'
        transcript_prediction_cuff(config_file)
    elif options.stringtie_pred:
        print 'Operation selected: Transcript assembly based on mapped RNASeq read data \
            with StringTie'
        transcript_prediction_stringtie(config_file)
    elif options.filter_trsk:
        print 'Operation selected: Filter out gene models from TranscriptSkimmer \
            predictions - criteria: splice-site consensus and length of the ORF.'
        filter_genes(config_file, "trsk")
    elif options.filter_cuff:
        print 'Operation selected: Filter out gene models from cufflinks predictions - \
            criteria: splice-site consensus, length of the ORF and read coverage to the \
            region.'
        filter_genes(config_file, "cufflinks")
    elif options.filter_db:
        print 'Operation selected: Filter out gene models from public database - criteria: \
            splice-site consensus and length of the ORF'
        filter_genes(config_file, "onlinedb")
    elif options.trsk_label:
        print 'Operation selected: Extract different genomic signal label sequences from \
            TranscriptSkimmer.'
        fetch_db_signals(config_file, "trsk")
    elif options.cuff_label:
        print 'Operation selected: Extract different genomic signal label sequences from \
            cufflinks.'
        fetch_db_signals(config_file, "cufflinks")
    elif options.db_label:
        print 'Operation selected: Extract different genomic signal label sequences from \
            online database files.'
        fetch_db_signals(config_file, "onlinedb")


def call_fetch_db_signals(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from signal_labels import generate_genome_seq_labels as fetch_labels 
    fasta_file, gff_file, signal_type, count, poslabels_cnt, neglabels_cnt, flank_nts, out_dir = args_list
    os.chdir(out_dir)
    fetch_labels.main(fasta_file, gff_file, signal_type, count, poslabels_cnt, neglabels_cnt, flank_nts)
    return "done" 


def fetch_db_signals(yaml_config, data_method):
    """
    get the genomic signal labels bases on the annotation from external database
    """
    operation_seleted = "6"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():
        if data_method == "trsk":
            gff_file = "%s/%s_trsk_genes.gff" % (det['read_assembly_dir'], org_name)
            out_dir = "%s/trsk_4K_labels" % det['labels_dir']## new label sequence dir 
        elif data_method == "cufflinks":
            gff_file = "%s/%s_cufflinks_genes.gff" % (det['read_assembly_dir'], org_name)
            out_dir = "%s/cuff_4K_labels" % det['labels_dir']
        elif data_method == "onlinedb":
            gff_file = "%s/%s_%s.gff" % (det['read_assembly_dir'], org_name, det['genome_release_db']) ## db_anno 
            out_dir = "%s/jmlr_1K_sm_labels" % det['labels_dir']
        
        if not os.path.isfile(gff_file):## check the file present or not  
            print "error: genome annotation file missing %s" % gff_file
            sys.exit(0)
       
        if not os.path.exists(out_dir): ## create the new label sequence dir 
            os.makedirs(out_dir)

        for the_file in os.listdir(out_dir): ## cleaning the existing one 
            file_path = os.path.join(out_dir, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception, e:
                print e 
    
        #import subprocess 
        ## get the label count for each organisms, essentially the max number of genes available 
        #cmd = "grep -P \"\tgene\t\" %s | wc -l" % gff_file
        #proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #count, err = proc.communicate() 
        #count = int(count.strip())
        
        ## depends on the genomic signal type 
        count = 5000
        signal_type = "tss"
        poslabels_cnt = 1000
        neglabels_cnt = 3000
        flank_nts = 1200 

        ## arguments to pygrid 
        arg = [[det['fasta'], gff_file, signal_type, count, poslabels_cnt, neglabels_cnt, flank_nts, out_dir]]
        job = pg.cBioJob(call_fetch_db_signals, arg) 

        ## native specifications 
        job.mem="5gb"
        job.vmem="5gb"
        job.pmem="5gb"
        job.pvmem="5gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "1:00:00"

        Jobs.append(job)
    print 
    print "sending genomic signal fetch jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


def call_filter_genes(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from rnaseq_align_assembly import refine_transcript_models as filter_tool 
    gtf_file, fasta_file, result_file = args_list
    filter_tool.filter_gene_models(gtf_file, fasta_file, result_file)
    return "done"


def filter_genes(yaml_config, data_method):
    """
    filter out invalid gene models from the provided genome annotation
    """
    operation_seleted = "f"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():
        if data_method == "cufflinks":
            gff_file = "%s/transcripts.gtf" % det['read_assembly_dir'] ## cufflinks run output file 
            outFile = "%s/%s_cufflinks_genes.gff" % (det['read_assembly_dir'], org_name) ## example: A_thaliana_cufflinks_genes.gff 
        elif data_method == "trsk":
            gff_file = "%s/tmp_trsk_genes.gff" % det['read_assembly_dir'] ## trsk run output file 
            outFile = "%s/%s_trsk_genes.gff" % (det['read_assembly_dir'], org_name) ## example: A_thaliana_trsk_genes.gff  
        else:
            gff_file = det['gtf'] ## public database genome annotation file 
            outFile = "%s/%s_%s.gff" % (det['read_assembly_dir'], org_name, det['genome_release_db']) ## example: A_thaliana_arabidopsis-tair10.gff  

        ## arguments to pygrid 
        arg = [[gff_file, det['fasta'], outFile]]
        job = pg.cBioJob(call_filter_genes, arg) 

        ## native specifications 
        job.mem="6gb"
        job.vmem="6gb"
        job.pmem="6gb"
        job.pvmem="6gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "2:00:00"

        Jobs.append(job)
    print 
    print "sending filter gene models jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


def call_transcript_prediction_cuff(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from rnaseq_align_assembly import transcript_assembly as trassembly
    org_db, num_threads = args_list
    trassembly.run_cufflinks(org_db, num_threads)
    return "done"


def transcript_prediction_cuff(yaml_config):
    """
    transcript prediction using cufflinks
    """
    operation_seleted = "c"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        arg = [[det, 4]]

        job = pg.cBioJob(call_transcript_prediction_cuff, arg) 

        ## native specifications 
        job.mem="96gb"
        job.vmem="96gb"
        job.pmem="24gb"
        job.pvmem="24gb"
        job.nodes = 1
        job.ppn = 4
        job.walltime = "32:00:00"

        Jobs.append(job)
    print 
    print "sending transcript assembly cufflinks jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


def call_transcript_prediction_trsk(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from rnaseq_align_assembly import transcript_assembly as trassembly
    org_db = args_list
    trassembly.run_trsk(org_db)
    return "done"

    
def transcript_prediction_trsk(yaml_config):
    """
    transcript prediction using TranscriptSkimmer
    """
    operation_seleted = "4"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        arg = [det]

        job = pg.cBioJob(call_transcript_prediction_trsk, arg) 

        ## native specifications 
        job.mem="32gb"
        job.vmem="32gb"
        job.pmem="32gb"
        job.pvmem="32gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "9:00:00"

        Jobs.append(job)
    print 
    print "sending transcript assembly trsk jobs to worker"
    print 

    local = False  ## cluster compute switch 
    processedJobs = pg.process_jobs(Jobs, local=local)


def call_transcript_prediction_stringtie(args_list):
    """
    wrapper for submitting jobs to pygrid
    """

    from rnaseq_align_assembly import transcript_assembly as tsa
    org_db = args_list
    tsa.run_stringtie(org_db)
    return "done"


def transcript_prediction_stringtie(yaml_config):
    """
    transcript prediction using StringTie
    """

    operation_seleted = "4"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        arg = [det]

        job = pg.cBioJob(call_transcript_prediction_stringtie, arg) 
        
        cpus = 1 
        ## native specifications 
        job.mem="12gb"
        job.vmem="12gb"
        job.pmem="12gb"
        job.pvmem="12gb"
        job.nodes = 1
        job.ppn = cpus
        job.walltime = "6:00:00"

        Jobs.append(job)
    print("\nsending transcript assembly stringtie jobs to worker\n")

    local_compute = True ## switching between local multithreading and cluster computing
    
    processedJobs = pg.process_jobs(Jobs, local=local_compute)


def find_uniq_reads(yaml_config):
    """
    find uniquely mapped reads from a bam file 
    """
    operation_seleted = "u"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)
    print "NOT YET IMPLEMENTED."
    sys.exit(0)


def call_alignment_filter(args_list):
    """
    wrapper for submitting jobs to pygrid
    """

    from rnaseq_align_assembly import star_align_rna as filter
    org_name, out_dir, num_cpus = args_list
    filter.run_mmr(org_name, out_dir, num_cpus)
    return "done"


def alignment_filter(yaml_config):
    """
    run multimapper resolution program 
    """
    operation_seleted = "m"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():

        num_cpus = 4
        ## arguments to pygrid 
        arg = [[det['short_name'], det['read_map_dir'], num_cpus]]

        job = pg.cBioJob(call_alignment_filter, arg) 

        ## native specifications 
        job.pmem="30gb"
        job.pvmem="30gb"
        job.mem="120gb"
        job.vmem="120gb"
        job.nodes = 1
        job.ppn = num_cpus
        job.walltime = "48:00:00"

        Jobs.append(job)
    print 
    print "sending multi map resolution jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


def call_align_reads(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from rnaseq_align_assembly import star_align_rna as rnastar 
    org_db, read_type, max_mates_gap_length, num_cpus = args_list
    rnastar.run_star_alignment(org_db, read_type, max_mates_gap_length, num_cpus) 
    return 'done'


def align_rnaseq_reads(yaml_config):
    """
    wrapper for aligning rnaseq reads using 
    """
    operation_seleted = "3"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        lib_type = 'PE'
        lib_type = 'SE' if len(det['fastq'])==1 else lib_type

        ## library insert size 
        lib_insert_size = 100000
        num_cpu = 4

        arg = [[det, lib_type, lib_insert_size, num_cpu]]

        job = pg.cBioJob(call_align_reads, arg) 
    
        job.mem="150gb"
        job.vmem="150gb"
        job.pmem="30gb"
        job.pvmem="30gb"
        job.nodes = 1
        job.ppn = num_cpu
        job.walltime = "48:00:00"
        
        Jobs.append(job)
    print 
    print "sending read alignment with STAR jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs, local=True)


def calculate_insert_size(yaml_config):
    """
    wrapper for calling calculate insert size function
    """
    operation_seleted = "i"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)
    print "NOT YET IMPLEMENTED."
    sys.exit(0)


def call_genome_index(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import prepare_data as ppd
    fasta_file, out_dir, genome_anno, num_workers, onematelength = args_list
    ppd.create_star_genome_index(fasta_file, out_dir, genome_anno, num_workers, onematelength)
    return 'done'


def create_genome_index(yaml_config):
    """
    wrapper for calling genome index function 
    """
    operation_seleted = "2"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        num_cpus = 4 
        arg = [[det['fasta'], det['genome_index_dir'], det['gtf'], num_cpus, det['read_length']-1]]

        job = pg.cBioJob(call_genome_index, arg) 
    
        job.mem="46gb"
        job.vmem="46gb"
        job.pmem="46gb"
        job.pvmem="46gb"
        job.nodes = 1
        job.ppn = num_cpus
        job.walltime = "24:00:00"
        
        Jobs.append(job)
    print 
    print "sending star genome index jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


def call_download_sra_file(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import download_data as dld
    sra_run_id, out_dir = args_list
    dld.download_sra_file(sra_run_id, out_dir)
    return 'done'


def download_sra_data(yaml_config):
    """
    download sra file for the working organism   
    """
    operation_seleted = "1"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = [] 
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        arg = [[det['sra_run_id'], det['fastq_path']]]

        job = pg.cBioJob(call_download_sra_file, arg) 
    
        job.mem="2gb"
        job.vmem="2gb"
        job.pmem="2gb"
        job.pvmem="2gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "1:00:00"
        
        Jobs.append(job)
    print 
    print "sending download SRA file jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


def call_decompose_sra_file(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import download_data as dld
    sra_file, out_dir = args_list
    dld.decompress_sra_file(sra_file, out_dir)
    return 'done'


def decompose_sra_file(yaml_config):
    """
    decompress the .sra file from ncbi sra
    """
    operation_seleted = "d"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = [] 
    for org_name, det in orgdb.items():
        sra_file = "%s/%s.sra"  % (det['fastq_path'], det['sra_run_id'])

        if not os.path.isfile(sra_file):## check the file present or not  
            print "error: missing sequencing read file %s" % sra_file
            sys.exit(0)
        
        ## TODO can be consider to the yaml file options 
        #library_type = "pe"
        library_type = "pe"
        compress_format = "gzip"

        ## arguments to pygrid 
        arg = [[sra_file, det['fastq_path']]]

        job = pg.cBioJob(call_decompose_sra_file, arg) 
    
        job.mem="6gb"
        job.vmem="6gb"
        job.pmem="6gb"
        job.pvmem="6gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "24:00:00"
        
        Jobs.append(job)
    print 
    print "sending decompress SRA file jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)

def call_fungi_fasta(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import download_data as dld
    release_num, organism, genome_path = args_list
    dld.fetch_ensembl_fungi_fasta(release_num, organism, genome_path)
    return 'done'
      
def call_metazoa_fasta(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import download_data as dld
    release_num, organism, genome_path = args_list
    dld.fetch_ensembl_metazoa_fasta(release_num, organism, genome_path)
    return 'done'

def call_phytozome_fasta(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import download_data as dld
    release_num, organism, genome_path = args_list
    dld.fetch_phytozome_fasta(release_num, organism, genome_path)
    return 'done'

def call_ensembl_fasta(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import download_data as dld
    release_num, organism, genome_path = args_list
    dld.fetch_ensembl_fasta(release_num, organism, genome_path)
    return 'done'

def download_fasta(yaml_config):
    """
    download fasta file from remote data publishing services
    """
    operation_seleted = "g"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = [] 
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        arg = [[det['release_nb'], det['long_name'], det['genome_dir']]]

        if det['release_db'] == 'ensembl_metazoa_genome':
            job = pg.cBioJob(call_metazoa_fasta, arg) 
        elif det['release_db'] == 'phytozome_genome':
            job = pg.cBioJob(call_phytozome_fasta, arg) 
        elif det['release_db'] == 'ensembl_genome':
            job = pg.cBioJob(call_ensembl_fasta, arg) 
        elif det['release_db'] == 'ensembl_fungi_genome':
            job = pg.cBioJob(call_fungi_fasta, arg) 
        else:
            exit("error: download fasta plugin for %s not available, module works with ensembl_genome, ensembl_metazoa_genome and phytozome_genome servers." % det['release_db'])

        job.mem="2gb"
        job.vmem="2gb"
        job.pmem="2gb"
        job.pvmem="2gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "2:00:00"
        
        Jobs.append(job)
    print 
    print "sending fasta download job to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


def call_fungi_gtf(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import download_data as dld
    release_num, organism, genome_path = args_list
    dld.fetch_ensembl_fungi_gtf(release_num, organism, genome_path)
    return 'done'

def call_metazoa_gtf(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import download_data as dld
    release_num, organism, genome_path = args_list
    dld.fetch_ensembl_metazoa_gtf(release_num, organism, genome_path)
    return 'done'

def call_phytozome_gtf(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import download_data as dld
    release_num, organism, genome_path = args_list
    dld.fetch_phytozome_gff(release_num, organism, genome_path)
    return 'done'

def call_ensembl_gtf(args_list):
    """
    wrapper for submitting jobs to pygrid
    """
    from fetch_remote_data import download_data as dld
    release_num, organism, genome_path = args_list
    dld.fetch_ensembl_gtf(release_num, organism, genome_path)
    return 'done'

def download_gtf(yaml_config):
    """
    download gtf/gff file from remote data publishing services
    """
    operation_seleted = "a"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = [] 
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        arg = [[det['release_nb'], det['long_name'], det['genome_dir']]]

        if det['release_db'] == 'ensembl_metazoa_genome':
            job = pg.cBioJob(call_metazoa_gtf, arg) 
        elif det['release_db'] == 'phytozome_genome':
            job = pg.cBioJob(call_phytozome_gtf, arg) 
        elif det['release_db'] == 'ensembl_genome':
            job = pg.cBioJob(call_ensembl_gtf, arg) 
        elif det['release_db'] == 'ensembl_fungi_genome':
            job = pg.cBioJob(call_fungi_gtf, arg) 
        else:
            exit("error: download gtf plugin for %s not available, module works with ensembl_genome, ensembl_metazoa_genome and phytozome_genome servers." % det['release_db'])

        job.mem="2gb"
        job.vmem="2gb"
        job.pmem="2gb"
        job.pvmem="2gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "2:00:00"
        
        Jobs.append(job)
    print 
    print "sending gtf download job to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


if __name__=="__main__":
    main() 
