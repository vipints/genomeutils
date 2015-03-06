#!/usr/bin/env python 
"""
master script to execute the pipeline
"""

import os 
import sys 
import yaml 

import subprocess 

import libpyjobrunner as pg

from optparse import OptionParser

from fetch_remote_data import download_data as dld

from signal_labels import experiment_details_db as expdb

assert sys.version_info[:2] >= ( 2, 4 )

def main():
    """
    Managing the experiment run in different levels 

    Options
    
    TODO    
    -1 download_public_data from SRA ENSEMBL Phytozome
        uncompressing the SRA file 
        manual cleaning of downloaded genome 
        create the genome genome indices
        calculate the insert size 

    -2 rnaseq read mapping process     
        multiple mapper issue 
        uniquely mapped reads 
    -3 transcript assembly  
        cufflinks prediction 
        filter gene models 
    -4 extract signal labels 
    """

    parser = OptionParser() 

    parser.add_option( "-1", "--download_public_data", action="store_true", dest="download_public_data", default=False, help="download public datasets" )
    parser.add_option( "-a", "--genome_index", action="store_true", dest="genome_index", default=False, help="Create STAR genome index to align the reads." )
    parser.add_option( "-2", "--read_mapping", action="store_true", dest="read_mapping", default=False, help="RNASeq read mapping to the genome using STAR." )
    parser.add_option( "-m", "--multi_map_resolve", action="store_true", dest="multi_map_resolve", default=False, help="Multimapper resolution (mmr) program on aligned reads." )
    parser.add_option( "-3", "--trsk_prediction", action="store_true", dest="trsk_prediction", default=False, help="Transcript assembly using TranscriptSkimmer." )
    parser.add_option( "-c", "--cufflinks_prediction", action="store_true", dest="cufflinks_prediction", default=False, help="Transcript assembly using Cufflinks." )
    parser.add_option( "-f", "--filter_trsk_out", action="store_true", dest="filter_trsk_out", default=False, help="Apply filter to the TRSK predicted gene models." )
    parser.add_option( "--filter_cuff_out", action="store_true", dest="filter_cuff_out", default=False, help="Apply filter to the cufflinks predicted gene models." )
    parser.add_option( "--filter_db_anno", action="store_true", dest="filter_db_anno", default=False, help="Apply filter to the online db gene models." )
    parser.add_option( "-4", "--fetch_trsk_labels", action="store_true", dest="fetch_trsk_labels", default=False, help="Fetch labels from TranscriptSkimmer." )
    parser.add_option( "--fetch_cuff_labels", action="store_true", dest="fetch_cuff_labels", default=False, help="Fetch labels from cufflinks." )
    parser.add_option( "--fetch_db_labels", action="store_true", dest="fetch_db_labels", default=False, help="Fetch labels from public database annotation files." )

    ( options, args ) = parser.parse_args()
    try:
        config_file = args[0]
    except:
        print __doc__
        sys.exit(-1)

    if not (options.download_public_data ^ options.genome_index ^ \
            options.read_mapping ^ options.multi_map_resolve ^ \
            options.trsk_prediction ^ options.cufflinks_prediction ^ \
            options.filter_trsk_out ^ options.fetch_trsk_labels ^ \
            options.filter_cuff_out ^ options.filter_db_anno ^ \
            options.fetch_cuff_labels ^ options.fetch_db_labels):
        parser.print_help()
        sys.exit(-1)
        
    print 'Using config file %s for the experiment.' % config_file

    if options.download_public_data:
        print 'Operation selected: Download public genome dataset and Sequencing experiment files'
        download_public_data(config_file)

    elif options.genome_index:
        print 'Operation selected: Create STAR genome index'
        create_genome_index(config_file)

    elif options.read_mapping: 
        print 'Operation selected: Read alignment with STAR'
        align_rnaseq_reads(config_file)

    elif options.multi_map_resolve:
        print 'Operation selected: Multiple read mapper resolution with MMR'
        alignment_filter(config_file) 

    elif options.trsk_prediction:
        print 'Operation selected: Transcript assembly based on mapped RNASeq read data with TranscriptSkimmer'
        transcript_prediction_trsk(config_file)

    elif options.cufflinks_prediction:
        print 'Operation selected: Transcript assembly based on mapped RNASeq read data with Cufflinks'
        transcript_prediction_cuff(config_file)

    elif options.filter_trsk_out:
        print 'Operation selected: Filter out gene models from TranscriptSkimmer predictions - criteria: splice-site consensus, length of the ORF and read coverage to the region.'
        filter_genes(config_file, "trsk")

    elif options.filter_cuff_out:
        print 'Operation selected: Filter out gene models from cufflinks predictions - criteria: splice-site consensus, length of the ORF and read coverage to the region.'
        filter_genes(config_file, "cufflinks")

    elif options.filter_db_anno:
        print 'Operation selected: Filter out gene models from public database - criteria: splice-site consensus, length of the ORF and read coverage to the region.'
        filter_genes(config_file, "onlinedb")

    elif options.fetch_trsk_labels:
        print 'Operation selected: Extract different genomic signal label sequences from TranscriptSkimmer.'
        fetch_db_signals(config_file, "trsk")

    elif options.fetch_cuff_labels:
        print 'Operation selected: Extract different genomic signal label sequences from cufflinks.'
        fetch_db_signals(config_file, "cufflinks")

    elif options.fetch_db_labels:
        print 'Operation selected: Extract different genomic signal label sequences from online database files.'
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

    operation_seleted = "4"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():

        if data_method == "trsk":
            gff_file = "%s/%s_trsk_genes.gff" % (det['read_assembly_dir'], org_name)
            out_dir = "%s/trsk_labels" % det['labels_dir']## new label sequence dir 
        elif data_method == "cufflinks":
            gff_file = "%s/%s_cufflinks_genes.gff" % (det['read_assembly_dir'], org_name)
            out_dir = "%s/cuff_labels" % det['labels_dir']
        else:
            gff_file = "%s/%s_%s.gff" % (det['read_assembly_dir'], org_name, det['genome_release_db']) ## db_anno 
            out_dir = "%s/db_labels" % det['labels_dir']
        
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
    
        ## get the label count for each organisms, essentially the max number of genes available 
        #cmd = "grep -P \"\tgene\t\" %s | wc -l" % gff_file
        #proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #count, err = proc.communicate() 
        #count = int(count.strip())

        count = 4000
        signal_type = "tss"
        poslabels_cnt = 1000
        neglabels_cnt = 3000
        flank_nts = 1200 

        ## arguments to pygrid 
        arg = [[det['fasta'], gff_file, signal_type, count, poslabels_cnt, neglabels_cnt, flank_nts, out_dir]]
        
        job = pg.cBioJob(call_fetch_db_signals, arg) 

        ## native specifications 
        job.mem="6gb"
        job.vmem="6gb"
        job.pmem="6gb"
        job.pvmem="6gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "3:00:00"

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
        arg = [[det, 18]]

        job = pg.cBioJob(call_transcript_prediction_cuff, arg) 

        ## native specifications 
        job.mem="72gb"
        job.vmem="72gb"
        job.pmem="18gb"
        job.pvmem="18gb"
        job.nodes = 1
        job.ppn = 4
        job.walltime = "18:00:00"

        Jobs.append(job)

    print 
    print "sending transcript assembly jobs to worker"
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
    operation_seleted = "3"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)
    Jobs = []
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        arg = [det]

        job = pg.cBioJob(call_transcript_prediction_trsk, arg) 

        ## native specifications 
        job.mem="16gb"
        job.vmem="16gb"
        job.pmem="16gb"
        job.pvmem="16gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "6:00:00"

        Jobs.append(job)

    print 
    print "sending transcript assembly jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


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
        ## arguments to pygrid 
        arg = [[det['short_name'], det['read_map_dir'], 12]]

        job = pg.cBioJob(call_alignment_filter, arg) 

        ## native specifications 
        job.mem="36gb"
        job.vmem="36gb"
        job.pmem="12gb"
        job.pvmem="12gb"
        job.nodes = 1
        job.ppn = 3
        job.walltime = "24:00:00"

        Jobs.append(job)

    print 
    print "sending jobs to worker"
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

    operation_seleted = "2"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        lib_type = 'PE'
        lib_type = 'SE' if len(det['fastq'])==1 else lib_type

        arg = [[det, lib_type, 100000, 16]]

        job = pg.cBioJob(call_align_reads, arg) 
    
        job.mem="48gb"
        job.vmem="48gb"
        job.pmem="12gb"
        job.pvmem="12gb"
        job.nodes = 1
        job.ppn = 4
        job.walltime = "24:00:00"
        
        Jobs.append(job)

    print 
    print "sending jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


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
    operation_seleted = "a"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)

    Jobs = []
    for org_name, det in orgdb.items():
        ## arguments to pygrid 
        arg = [[det['fasta'], det['genome_index_dir'], det['gtf'], 12, det['read_length']-1]]

        job = pg.cBioJob(call_genome_index, arg) 
    
        job.mem="32gb"
        job.vmem="32gb"
        job.pmem="32gb"
        job.pvmem="32gb"
        job.nodes = 1
        job.ppn = 1
        job.walltime = "16:00:00"
        
        Jobs.append(job)

    print 
    print "sending jobs to worker"
    print 
    processedJobs = pg.process_jobs(Jobs)


def download_public_data(yaml_config):
    """
    """
    operation_seleted = "1"
    orgdb = expdb.experiment_db(yaml_config, operation_seleted)
    print orgdb


    #import ipdb 
    #ipdb.set_trace()



def download_fasta(org_details):
    """
    download fasta file from remote data publishing services
    """

    for org_name, det in org_details.items():
        if det['release_db'] == 'ensembl':
            dld.fetch_ensembl_fasta(det['release_num'], det['name'], det['fasta'])
        elif det['release_db'] == 'phytozome':
            dld.fetch_phytozome_fasta(det['release_num'], det['name'], det['fasta']) 
        elif det['release_db'] == 'ensembl_metazoa':
            dld.fetch_ensembl_metazoa_fasta(det['release_num'], det['name'], det['fasta'])
        else:
            print "download fasta plugin for %s not available, module works with ensembl, ensembl_metazoa and phytozome." % det['release_db']


def download_gtf(org_details):
    """
    download gtf/gff file from remote data publishing services
    """
    for org_name, det in org_details.items():
        if det['release_db'] == 'ensembl':
            dld.fetch_ensembl_gtf(det['release_num'], det['name'], det['gtf'])
        elif det['release_db'] == 'phytozome':
            dld.fetch_phytozome_gff(det['release_num'], det['name'], det['gtf']) 
        elif det['release_db'] == 'ensembl_metazoa':
            dld.fetch_ensembl_metazoa_gtf(det['release_num'], det['name'], det['gtf'])
        else:
            print "download gtf plugin for %s not available, module works with ensembl, ensembl_metazoa and phytozome." % det['release_db']


def download_uncompress_sra_file(org_details):
    """
    download and uncompress the sra files 
    """

    for org_name, det in org_details.items():
        sra_file = dld.download_sra_file(det['sra_run_id'], det['fastq_path'])
        
        library_type = "pe"
        compress_format = "gzip"
        dld.uncompress_sra_file(sra_file, det['fastq_path'], library_type, compress_format)


if __name__=="__main__":
    
    main() 

    """
    from signal_labels import org_details_db as odb
    org_details = odb.make_org_db(infile, data_path, exp_path) 

    download_fasta(org_details) 

    download_gtf(org_details)     

    download_uncompress_sra_file(org_details)

    #TODO
    save a pickle file with organisms details with updated genome annotation and sra file, path informations 

    the object will be passed to the next preprocessing manual tweeking

    
    #FIXME 
    fas_file = ""
    chr_names = prd.read_genome_file(fas_file) 

    fas_out = "" 
    prd.clean_genome_file(chr_names, fas_file, fas_out):

    gtf_out = "" 
    prd.clean_anno_file(chr_names, gtf_file, gtf_out)
    """
