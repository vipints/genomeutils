#!/usr/bin/env python 
"""
wrapper program to run transcriptome assembly using TransriptSkimmer 
and cufflinks on sequencing read alignment data. 

Requirement: 
    TransriptSkimmer - 
    cufflinks - 

    python libraries: 
    pysam 
"""

import os 
import sys
import pysam 
import shutil
import subprocess


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
            
        bam_file = "%s.bam" % sorted_bam

    print 'using bam file from %s' % bam_file
    if not os.path.exists(bam_file + ".bai"):
        sys.stdout.write('indexing... ')
        pysam.index(bam_file) 
        sys.stdout.write(' ...done\n')

    ## always use quiet mode to avoid problems with storing log output.
    cli_cuff = "cufflinks -q --no-update-check \
        -F %.2f \
        -I %d \
        --min-intron-length %d \
        --library-type fr-unstranded \
        -p %d \
        -o %s \
        %s" % (min_isoform_frac, max_intron_length, min_intron_length, num_cpus, result_dir, bam_file)
  
    sys.stdout.write('\trun cufflinks as: %s \n' % cli_cuff)
    
    try:
        os.chdir(result_dir)
        ## run the command
        process = subprocess.Popen(cli_cuff, shell=True) 
        process.wait()

    except Exception, e:
        print 'Error running cufflinks.\n%s' %  str( e )
        

def run_trsk(org_db, out_gff_file="tmp_trsk_genes.gff"):
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
            
        bam_file = "%s.bam" % sorted_bam
        
    print 'using bam file from %s' % bam_file
    if not os.path.exists(bam_file + ".bai"):
        sys.stdout.write('indexing... ')
        pysam.index(bam_file) 
        sys.stdout.write(' ...done\n')

    max_intron_length = org_db['max_intron_len']
    max_exon_length = org_db['max_exon_len'] 
    ##FIXME to be included in the genome database  
    max_intergenic_region = 10000 

    result_dir = org_db['read_assembly_dir']
    gio_path_temp = os.path.join(result_dir, "temp_gio")
    make_gio(genome_seq_file, gio_path_temp)
    gio_file = "%s/genome.config" % gio_path_temp 

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

    cli_trsk = "infer_genes -gio %s -bam %s -gff %s %s" % (gio_file, bam_file, out_gff_file, options)  
    sys.stdout.write('\trun TransriptSkimmer as: %s \n' % cli_trsk)
    
    try:
        os.chdir(result_dir)
        process = subprocess.Popen(cli_trsk, shell=True) 
        process.wait()

    except Exception, e:
        print 'Error running TranscriptSkimmer.\n%s' %  str( e )
   
    ## cleaning 
    shutil.rmtree(gio_path_temp)


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
		

