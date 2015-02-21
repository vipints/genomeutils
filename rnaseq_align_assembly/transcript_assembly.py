#!/usr/bin/env python 
"""
Program to run transcriptome assembly using TransriptSkimmer and cufflinks program on read alignment data 

Requirement: 
    TransriptSkimmer - 
    cufflinks - 

Standard python libraries: 
    numpy 
    pysam 
    biopython 
    gfftools 

    TODO sync the validate_pred_gene_models function to the refine_ module 
"""

from __future__ import division
import os 
import sys
import shutil
import subprocess

import pysam 
import numpy 

import collections
from Bio import SeqIO 
from gfftools import GFFParser, helper 


def run_cufflinks(org_db, num_cpus=4):
    """
    run cufflinks program on mapped reads 
    """
     
    org_name = org_db['short_name'] 
    print "preparing for cufflinks run for organism %s" % org_name

    min_isoform_frac = 0.25
    max_intron_length = org_db['max_intron_len']
    min_intron_length = 20

    result_dir = org_db['read_assembly_dir']

    bam_file = "%s/%s_Aligned_mmr_sortbyCoord.bam" % (org_db['read_map_dir'], org_name)
    if not os.path.isfile(bam_file):
        print "failed to fetch sorted mmr BAM file for organism: %s, trying to get the mmr file..." % org_name
        
        bam_file = "%s/%s_Aligned_mmr.bam" % (org_db['read_map_dir'], org_name)
        if not os.path.isfile(bam_file):
            print "error: failed to fetch mmr BAM file for organism %s" % org_name
            sys.exit(-1)
        
        ## sorting, indexing the bam file 
        file_prefix, ext = os.path.splitext(bam_file)
        sorted_bam = "%s_sortbyCoord" % file_prefix

        print "trying to sort based by the coordinates with output prefix as: %s" % sorted_bam
        if not os.path.isfile("%s.bam" % sorted_bam):
            print 'sorting...'
            pysam.sort(bam_file, sorted_bam)
            
        sorted_bam = "%s.bam" % sorted_bam
        print "now creating the index for %s " % sorted_bam
        if not os.path.exists(sorted_bam + ".bai"):
            print 'indexing...'
            pysam.index(sorted_bam) 

        bam_file = sorted_bam 
        print "done"

    print "using bam file from %s" % bam_file

    ## always use quiet mode to avoid problems with storing log output.
    cli_cuff = "cufflinks -q --no-update-check \
        -F %.2f \
        -I %d \
        --min-intron-length %d \
        --library-type fr-unstranded \
        -p %d \
        -o %s \
        %s" % (min_isoform_frac, max_intron_length, min_intron_length, num_cpus, result_dir, bam_file)
  
    sys.stdout.write('\trun cufflinks as %s \n' % cli_cuff)
    
    try:
        os.chdir(result_dir)
        ## run the command
        process = subprocess.Popen(cli_cuff, shell=True) 
        process.wait()

        ## cleaning 
        print "cleaning the predicted transcript models"
        out_file = "%s_cufflinks_genes.gff" % org_name
        
        genome_seq_file = org_db['fasta']
        
        out_cuff_gtf = "%s/transcripts.gtf" % result_dir
        ## transcripts.gtf default file for cufflinks prediction 
        final_transcripts = validate_pred_gene_models(out_cuff_gtf, genome_seq_file, out_file)

        print 'cufflinks predicted transcripts are stored at %s/%s' % (result_dir, final_transcripts) 
        os.unlink("%s/%s" % (result_dir, "transcripts.gtf"))

    except Exception, e:
        print 'Error running cufflinks.\n%s' %  str( e )
        
    

def run_trsk(org_db, out_gff_file="_tmp_trsk_genes.gff"):
    """
    run TransriptSkimmer with mapped reads and genome sequence
    """

    org_name = org_db['short_name'] 
    print "preparing for TransriptSkimmer run for organism %s" % org_name

    genome_seq_file = org_db['fasta']

    if not os.path.isfile(genome_seq_file):
        print 'error: failed to fetch genome sequence file %s for organism %s' % (genome_seq_file, org_name) 
        sys.exit(-1)

    print "using genome sequence file from %s" % genome_seq_file

    ## expect the mmr result file in below format ex: A_thaliana/read_mapping/A_thaliana_Aligned_mmr_sortbyCoord.bam 
    bam_file = "%s/%s_Aligned_mmr_sortbyCoord.bam" % (org_db['read_map_dir'], org_name)
    if not os.path.isfile(bam_file):
        print "failed to fetch sorted mmr BAM file for organism: %s, trying to get the mmr file..." % org_name
        
        bam_file = "%s/%s_Aligned_mmr.bam" % (org_db['read_map_dir'], org_name)
        if not os.path.isfile(bam_file):
            print "error: failed to fetch mmr BAM file for organism %s" % org_name
            sys.exit(-1)
        
        ## sorting, indexing the bam file 
        file_prefix, ext = os.path.splitext(bam_file)
        sorted_bam = "%s_sortbyCoord" % file_prefix

        print "trying to sort based by the coordinates with output prefix as: %s" % sorted_bam
        if not os.path.isfile("%s.bam" % sorted_bam):
            print 'sorting...'
            pysam.sort(bam_file, sorted_bam)
            
        sorted_bam = "%s.bam" % sorted_bam
        print "now creating the index for %s " % sorted_bam
        if not os.path.exists(sorted_bam + ".bai"):
            print 'indexing...'
            pysam.index(sorted_bam) 

        bam_file = sorted_bam 
        print "done"

    print 'using bam file from %s' % bam_file

    max_intron_length = org_db['max_intron_len']
    max_exon_length = org_db['max_exon_len'] 
    ##FIXME to be included in the genome database  
    max_intergenic_region = 10000 

    result_dir = org_db['read_assembly_dir']
    gio_path_temp = os.path.join(result_dir, "temp_gio")
    make_gio(genome_seq_file, gio_path_temp)

    options="-maxel %d \
    -nss \
    -reglen 0.66 \
    -maxic %d \
    -minic 20 \
    -maxin %d \
    -mm 4 \
    -exm 3 \
    -indt 150 \
    -exd 20 \
    -tf 0.5 \
    -inscf 3 \
    -excut 3 \
    -toff 100 \
    -el 15" % (max_exon_length, max_intergenic_region, max_intron_length)

    os.chdir(result_dir)

    gio_file = "%s/genome.config" % gio_path_temp 

    cli_trsk = "infer_genes -gio %s -bam %s -gff %s %s" % (gio_file, bam_file, out_gff_file, options)  
    sys.stdout.write('\trun TransriptSkimmer as %s \n' % cli_trsk)

    process = subprocess.Popen(cli_trsk, shell=True) 
    process.wait()
    
    ## cleaning 
    print "cleaning the predicted transcript models"
    out_file = "%s_trsk_genes.gff" % org_name
    final_transcripts = validate_pred_gene_models(out_gff_file, genome_seq_file, out_file)
    print "predicted transcript models are available at: %s/%s" % (result_dir, final_transcripts)  

    shutil.rmtree(gio_path_temp)
    os.unlink("%s/%s" % (result_dir, out_gff_file))


def make_gio(in_file_name, gio_path):
	"""
    make_gio builds a genome information object for an input fasta file. 

	takes 2 arguments:

	@args fasta_file: is the input file in fasta format
    @type fasta_file: str 
	@args gio_path: is the directory to which the genome information object will be written to
    @type gio_path: dir 
	"""

	try:
		f_in = file(in_file_name, "r")
	except Exception, msg:
		print msg
		print "cannot open infile '" + in_file_name + "'"
		sys.exit(1)
	   
	write_dna = 1 
	flat_path = os.path.join(gio_path, "genome")
	try:
		if os.path.exists(flat_path):
			print "directory " + flat_path + " exists already."
		else:
			os.makedirs(flat_path)

	except Exception, msg:
		print msg
		print "cannot create path '" + flat_path + "'"
		sys.exit(1)

	f_out = None
	f_out_dna = None
	contig_list = []
	num_bases = 0
	
	for line in f_in:

		if line.isspace():
			print "warning: wrong format. ignoring empty line in file '" + in_file_name + "'"
			continue

		if line[0].isspace():
			print "wrong format: leading white space in file '" + in_file_name + "'"
			sys.exit(1)
	
		if line.startswith(">"):
			
			if f_out != None:
				f_out.close()

			if f_out_dna != None:
				f_out_dna.close()

			contig_list.append(line[1:-1].split()[0])
			out_name = os.path.join(flat_path, contig_list[-1] + ".flat")
			out_dna_name = os.path.join(flat_path, contig_list[-1] + ".dna")
			try:
				f_out = file(out_name, "w")
				#print "creating file '" + out_name + "'"
				if write_dna==1:
					f_out_dna = file(out_dna_name, "w")
					f_out_dna.write(line)
					#print "creating file '" + out_dna_name + "'"

			except Exception, msg:
				print msg
				print "cannot open file '" + out_name + "'"
				sys.exit(1)
				
		else:
			try:
				f_out.write(line[0:-1].lower())
				if write_dna==1:
					f_out_dna.write(line.lower())
			except Exception, msg:
				if f_out != None:
					print msg
					print "cannot write to file '" +out_name + "'"
					sys.exit(1)
				else:
					print "improper input format. No header in first line"
					sys.exit(1)

			num_bases += len(line)-1

	f_out.close()

	try:
		print "creating file '" + os.path.join(gio_path, "genome.config") + "'"
		f_conf = file(os.path.join(gio_path, "genome.config"),"w")
		f_conf.write("BASEDIR " +  os.path.abspath(gio_path) +"\n\n")
		f_conf.write("CONTIGS " +  str(len(contig_list)) +"\n")
		for c in contig_list:
			f_conf.write(c + "\tgenome/" + c + ".flat\tgenome/" + c + ".dna\n")
		f_conf.write("\nALPHABET acgt\n\n")
		f_conf.write("ESTFILES 0\n\n")
		f_conf.write("CDNAFILES 0\n\n")
		f_conf.write("ANNOTATIONFILES 0\n")
		f_conf.close()
	except Exception, msg:
		print msg
		print "cannot create file '" + os.path.join(gio_path, "genome.config") + "'"
		sys.exit(1)
		

def validate_pred_gene_models(gff_name, fas_file, out_fname="trsk_genes.gff"):
    """
    check the sequence consistency/quality of predicted fragment

    @args gff_name: result file gff format from TranscriptSkimmer
    @type gff_name: str
    @args fas_file: genome sequence in fasta format
    @type fas_file: str 
    @args out_fname: filtered gene models in gff format (default: trsk_genes.gff)
    @type out_fname: str 
    """

    ## getting the genome annotation from GFF file 
    gff_content = GFFParser.Parse(gff_name)
    
    ## getting the spliced transcripts from the predicted gene list 
    transcripts_region = collections.defaultdict(list)
    for gene_recd in gff_content:
        spliced_transcript = collections.defaultdict(list)

        for idx, sub_rec in enumerate(gene_recd['transcripts']):
            exon_cnt = len(gene_recd['exons'][idx])

            ## skipping the single-exon transcripts 
            if exon_cnt > 1: 
                for idk, ex in enumerate(gene_recd['exons'][idx]):
                    if idk == 0:
                        ex[0] = None 
                    if exon_cnt-1 == idk:
                        ex[1] = None

                    spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append(ex)

        transcripts_region[gene_recd['chr']].append(spliced_transcript)

    print "check for the splice site consensus for predicted transcripts"
    ## check for splice site consensus sequence of predicted transcripts 
    get_gene_models = collections.defaultdict()
    for fas_rec in SeqIO.parse(fas_file, "fasta"):
        if fas_rec.id in transcripts_region:
            for details in transcripts_region[fas_rec.id]:
                for genes, regions in details.items():

                    acc_cons_cnt = 0 
                    don_cons_cnt = 0 

                    for region in regions:
                        if genes[-1] == '+':
                            ## acceptor splice site 
                            if not numpy.isnan(region[0]):
                                acc_seq = fas_rec.seq[int(region[0])-3:int(region[0])-1]
                                if str(acc_seq).upper() == "AG":
                                    acc_cons_cnt += 1 

                            if not numpy.isnan(region[1]):
                                don_seq = fas_rec.seq[int(region[1]):int(region[1])+2]
                                if str(don_seq).upper() == "GT":
                                    don_cons_cnt +=1 

                        elif genes[-1] == '-':
                            ## donor splice site 
                            if not numpy.isnan(region[0]):
                                don_seq = fas_rec.seq[int(region[0])-3:int(region[0])-1]
                                don_seq = don_seq.reverse_complement()
                                if str(don_seq).upper() == "GT":
                                    don_cons_cnt +=1 
                            
                            if not numpy.isnan(region[1]):
                                acc_seq = fas_rec.seq[int(region[1]):int(region[1])+2]
                                acc_seq = acc_seq.reverse_complement()
                                if str(acc_seq).upper() == "AG":
                                    acc_cons_cnt += 1 
                    ## check for half of the consensus sites 
                    if acc_cons_cnt > (len(regions)/2) and don_cons_cnt > (len(regions)/2):
                        get_gene_models[(fas_rec.id, genes[0], genes[1], genes[2])] = 1   
    
    gff_cont = GFFParser.Parse(gff_name)

    ## filter out the best gene models based on the consensus 
    print "writing the fine tuned transctipts to the the file"
    out_fh = open(out_fname, "w")
    for recd in gff_cont:
        trans_indices = [] 

        for idx, sub_rec in enumerate(recd['transcripts']):
            if (recd['chr'], recd['name'], sub_rec[0], recd['strand']) in get_gene_models:
                trans_indices.append(idx)

        if trans_indices:
            chr_name = recd['chr']
            strand = recd['strand']
            start = recd['start']
            stop = recd['stop']
            source = recd['source']
            ID = recd['name']
            Name = recd['gene_info']['Name']
            Name = ID if Name != None else Name  
            out_fh.write('%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s\n' % (chr_name, source, start, stop, strand, ID, Name))
                
            for idz, tid in enumerate(recd['transcripts']):
                if idz in trans_indices:

                    t_start = recd['exons'][idz][0][0]
                    t_stop = recd['exons'][idz][-1][-1]
                    t_type = recd['transcript_type'][idz] 

                    out_fh.write('%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (chr_name, source, t_type, t_start, t_stop, strand, tid[0], ID))
                    
                    for ex_cod in recd['utr5_exons'][idz]:
                        out_fh.write('%s\t%s\tfive_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s\n' % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0])) 
                    for ex_cod in recd['cds_exons'][idz]:
                        out_fh.write('%s\t%s\tCDS\t%d\t%d\t.\t%s\t%d\tParent=%s\n' % (chr_name, source, ex_cod[0], ex_cod[1], strand, ex_cod[2], tid[0])) 
                    for ex_cod in recd['utr3_exons'][idz]:
                        out_fh.write('%s\t%s\tthree_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s\n' % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0]))
                    for ex_cod in recd['exons'][idz]:
                        out_fh.write('%s\t%s\texon\t%d\t%d\t.\t%s\t.\tParent=%s\n' % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0])) 
    out_fh.close()
    return out_fname
