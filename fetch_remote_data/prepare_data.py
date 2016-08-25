#!/usr/bin/env python 
"""
modules for large scale experiment runs includes:
    - reading of genome sequence in fasta format 
    - manual cleaning of genome and annotation files
    - feature annotation db for querying details
    - creating star genome indicies

Requirement:
    STAR aligner: https://github.com/alexdobin/STAR 
    gfftools: https://github.com/vipints/genomeutils/tree/master/gfftools 
    Biopython: http://biopython.org
"""

import os 
import re
import sys 
import shutil
import subprocess 
from Bio import SeqIO 
from gfftools import helper, GFFParser


def stop_err(msg):
    """
    stop the execution and print out the captured error message. 
    """
    sys.stderr.write('%s\n' % msg)
    sys.exit(-1)


def clean_anno_file(chr_names, gtf_file, gtf_out):
    """
    make stable annotation file with valid contig name 

    @args chr_names: different contig names with a valid genome sequence 
    @type chr_names: dict
    @args gtf_file: genome annotation in gtf/gff form 
    @type gtf_file: str 
    @args gtf_out: new genome annotation in gtf/gff form 
    @type gtf_out: str 
    """

    # get the filehandler from input file
    try:
        fh = helper.open_file(gtf_file)
    except Exception, errmsg:
        stop_err('error %s in reading file %s' % (errmsg, gtf_file)) 

    # check the out filehandler
    try:
        outfh = open(gtf_out, "w")
    except Exception, errmsg:
        stop_err('error %s in writing file %s' % (errmsg, gtf_out)) 

    for line in fh:
        line = line.strip('\n\r')

        ## preserving the fasta header if present 
        if line[0] in  ['#', '>']:
            outfh.write(line + '\n')
            continue
        ## preserving the genome sequence if present
        if not re.search('\t', line):
            outfh.write(line + '\n')
            continue
        ## looking for gtf/gff files 
        fields = line.split('\t')
        assert len(fields) >= 8, fields
         
        if fields[0] in chr_names:
            outfh.write(line + '\n')

    fh.close() 
    outfh.close()


def read_genome_file(fas_file):
    """
    read genome file in fasta and return the list of chromosomes/contigs 

    @args fas_file: genome sequence in fasta file  
    @type fas_file: str 

    returns a list with contig_names and length
    """
    
    # get the filehandler from input file
    try:
        fh = helper.open_file(fas_file)
    except Exception, errmsg:
        stop_err('error in reading file '+ errmsg) 

    chrom_names = [] 
    for rec in SeqIO.parse(fh, "fasta"):
        print "parsing contig %s details" % rec.id 
        chrom_names.append((rec.id, len(rec.seq)))
    
    fh.close()
    # return the list with chromosome identifier and its sequence length 
    return chrom_names

    """
    based on eye inspection, the returned list can be trimmed and 
    create a dictionary with the best chromosomes , something like: 

    Take the list 0-15
    chr_best = chrom_names[0:15]
    
    change to dict
    chr_best = dict(chr_best) 

    and finally this dictionary can be passed to the genome cleaning 
    function - clean_genome_file
    """


def clean_genome_file(chr_names, fas_file, fas_out):
    """
    make a stable genome file with valid contigs 

    @args chr_names: different contig names with a valid genome sequence 
    @type chr_names: dict
    @args fas_file: genome sequence in fasta file  
    @type fas_file: str 
    @args fas_out: new genome sequence file in fasta format 
    @type fas_out: str 
    """

    # get the filehandler from input file
    try:
        fh = helper.open_file(fas_file)
    except Exception, errmsg:
        stop_err('error in reading file '+ errmsg) 

    # check the out filehandler
    try:
        outfh = open(fas_out, "w")
    except Exception, errmsg:
        stop_err('error in writing file '+ errmsg) 

    # writing stable contig genome sequence in FASTA format 
    for rec in SeqIO.parse(fh, "fasta"):
        if rec.id in chr_names:
            outfh.write(rec.format("fasta"))
            print "writing the contig %s details" % rec.id  

    fh.close()
    outfh.close()


def genome_file_rec_extract(chr_pattn, fas_file, fas_out):
    """
    get all contings based on a matiching string in the record identifier

    @args chr_pattn: pattern to be searched in contig names 
    @type chr_pattn: str
    @args fas_file: genome sequence in fasta file  
    @type fas_file: str 
    @args fas_out: new genome sequence file in fasta format 
    @type fas_out: str 
    """

    # get the filehandler from input file
    try:
        fh = helper.open_file(fas_file)
    except Exception, errmsg:
        stop_err('error in reading file '+ errmsg) 

    # check the out filehandler
    try:
        outfh = open(fas_out, "w")
    except Exception, errmsg:
        stop_err('error in writing file '+ errmsg) 

    # writing stable contig genome sequence in FASTA format 
    for rec in SeqIO.parse(fh, "fasta"):
        if re.search(chr_pattn, rec.id):
            outfh.write(rec.format("fasta"))
            print "writing the contig %s details" % rec.id  

    fh.close()
    outfh.close()


def make_anno_db(gff_file): 
    """
    extract the features from a gtf/gff file and store efficiently to query 

    @args gff_file: genome annotation file
    @type gff_file: str 
    """

    gff_cont = GFFParser.Parse(gff_file)  

    intron_size = dict() 
    exon_size = dict() 

    for rec in gff_cont:
        for idx, tid in enumerate(rec['transcripts']):

            if not rec['exons'][idx].any():
                continue

            try: # (Pdb) rec['exons'][0] -> array(nan)
                import numpy as np 
                if np.isnan(rec['exons'][idx]):
                    continue
            except:
                pass 
                    
            try:
                exon_cnt = len(rec['exons'][idx])
            except:
                continue
                
            if exon_cnt > 1:
                intron_start = 0 
                
                for xq, excod in enumerate(rec['exons'][idx]): 
                    
                    if xq > 0: 
                        #print intron_start, excod[0]-1 
                        if excod[0]-intron_start==1:
                            intron_start = excod[1]+1
                            exon_size[intron_start-excod[0]] = 1
                            continue

                        intron_size[excod[0]-intron_start] = 1 
                        #print excod[0]-intron_start

                    intron_start = excod[1]+1
                    exon_size[intron_start-excod[0]] = 1
                    #print intron_start-excod[0]
    
    feat_db = dict() 
    if intron_size:
        keys_int = sorted(intron_size)
        keys_ex = sorted(exon_size)
        #print 'MaxIntronLength %d %d %d'  %(keys_int[-1], keys_int[-2], keys_int[-3])
        feat_db['min_intron'] = int(keys_int[0])
        feat_db['max_intron'] = int(keys_int[-3])

        feat_db['min_exon'] = int(keys_ex[0])
        feat_db['max_exon'] = int(keys_ex[-3])
        #print 'MaxExonLength %d %d %d'  %(keys_ex[-1], keys_ex[-2], keys_ex[-3]) 

        return feat_db 
    else:
        print "Error in feature mapping in file %s, please check the source of parent child features" % gff_file
        sys.exit(-1)


def create_star_genome_index(fasta_file, out_dir, genome_anno=None, num_workers=1, onematelength=100):
    """
    Creating STAR genome index with or without using genome annotation

    @args fasta_file: reference genome sequence file .fasta format 
    @type fasta_file: str 
    @args out_dir: genome index binary file storage place  
    @type out_dir: str 
    @args genome_anno: genome annotation file (optional) 
    @type genome_anno: str 
    @args num_workers: number of threads to run (default value = 1)
    @type num_workers: int 
    @args onematelength: One Mate Length (default value=100) 
    @type onematelength: int 
    """

    try:
        subprocess.call(["STAR"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        exit("Please make sure that the `STAR` binary is in your $PATH")
    
    file_prefx, ext = os.path.splitext(fasta_file)
    if ext in [".bz2", ".gz", ".lzma"]: ## checking for the compressed form of the file extension 
        exit("error: STAR - Generating genome indexes - recommended to use the uncompressed FASTA file %s." % fasta_file)
    
    if not genome_anno:
        cli_cmd = 'STAR \
        --runMode genomeGenerate \
        --genomeDir %s \
        --genomeFastaFiles %s \
        --runThreadN %d' % (out_dir, fasta_file, num_workers) 
    else:
        file_prefx, ext = os.path.splitext(genome_anno)
        if ext in [".bz2", ".gz", ".lzma"]: 
            exit("error: STAR - Generating genome indexes - recommended to use the uncompressed GTF/GFF file %s." % genome_anno)

        ## check for the file type  
        gff_hand = helper.open_file(genome_anno)
        for rec in gff_hand:
            rec = rec.strip('\n\r')

            # skip empty line fasta identifier and commented line
            if not rec or rec[0] in  ['#', '>']:
                continue
            # skip the genome sequence 
            if not re.search('\t', rec):
                continue
            parts = rec.split('\t')
            assert len(parts) >= 8, rec

            ftype, tags = GFFParser.attribute_tags(parts[-1])
            break 
        gff_hand.close() 

        ## according to the file type 
        if ftype:
            cli_cmd = 'STAR \
            --runMode genomeGenerate \
            --genomeDir %s \
            --genomeFastaFiles %s \
            --runThreadN %d \
            --sjdbGTFfile %s \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbOverhang %d' % (out_dir, fasta_file, num_workers, genome_anno, onematelength) 
        else:
            cli_cmd = 'STAR \
            --runMode genomeGenerate \
            --genomeDir %s \
            --genomeFastaFiles %s \
            --runThreadN %d \
            --sjdbGTFfile %s \
            --sjdbGTFfeatureExon exon \
            --sjdbOverhang %d' % (out_dir, fasta_file, num_workers, genome_anno, onematelength) 

    ## create downloadpath if doesnot exists 
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except OSError:
            exit("error: cannot create the directory %s." % out_dir)
    else:## if present any other old index files clean up the folder 
        for the_file in os.listdir(out_dir):
            file_path = os.path.join(out_dir, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception, e:
                print(e) 
    
    ## start the indexing job 
    sys.stdout.write('\trunning program as: %s \n' % cli_cmd)
    try:
        ## changing the working dir to run STAR 
        os.chdir(out_dir)
        ## Run command.
        process = subprocess.Popen(cli_cmd, shell=True) 
        returncode = process.wait()

        ## Error checking.
        if returncode != 0:
            raise Exception, "return code = %i" % returncode
        
        print("\nGenome index files are stored at %s\n" % out_dir)

    except Exception, e:
        exit('Error running STAR.\n%s' %  str( e ))

def parse_list(line, nb_elts):
    """
    specific to BWA-MEM stderr format
    """
    return map(lambda x: int(float(x)), ' '.join(line.strip().replace(',','').split()[-nb_elts:])[1:-1].split())


def calculate_insert_size_fastq(ref_genome, fastq_dir, fq_files):
    """
    calculate the library insert size from raw read files and reference genome sequence

    @args ref_genome: genome sequence file 
    @type ref_genome: str
    @args fastq_dir: fastq directory for the experiment 
    @type fastq_dir: str 
    @args fq_files: fastq file names in python list
    @type fq_files: list 
    """
    
    try:
        subprocess.call(["bwa"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        exit("Please make sure that the `bwa` binary is in your $PATH")
    
    ## check for the bwa index  
    if not os.path.isfile("%s.sa" % ref_genome):
        try:
            cmd_idx = "bwa index %s" % ref_genome 
            sys.stdout.write("bwa indexing starting...\n")
            process = subprocess.Popen(cmd_idx, shell=True) 
            returncode = process.wait()

            if returncode !=0:
                raise Exception, "Exit status return code = %i" % returncode
            sys.stdout.write("bwa indexing finished\n")

        except Exception, e:
            exit('Error running bwa index.\n%s' %  str( e ))

    ##adapted from: 
    ##Quickly estimates insert sizes of read datasets, given some sequence(s) they can be mapped to.
    ##Author: Rayan Chikhi         
    ##example: 
    ##     estimate-insert-sizes contigs.fa readsA_1.fq readsA_2.fq readsB_1.fq readsB_2.fq
    if len(fq_files) > 1:
        nb_cpu = 1 
        zip_type = {".gz" : "gzip -c", ".bz2" : "bzip2 -d -c"} 
        file_prefx, ext = os.path.splitext(fq_files[0])

        stats = dict()
        read1 = "%s/%s" % (fastq_dir, fq_files[0])
        read2 = "%s/%s" % (fastq_dir, fq_files[1])
        print("processing : \n %s \n %s " % (read1, read2))

        cmd = ["bwa", "mem"] + ["-t %d" % nb_cpu, ref_genome] + ["<(%s %s)" % (zip_type[ext], read1)] \
            +["<(%s %s)" % ( zip_type[ext], read2)]
        DEVNULL = open(os.devnull, 'wb')
        process = subprocess.Popen(cmd, stdout=DEVNULL, stderr=subprocess.PIPE)
        seen_candidate_line = False

        while True:
            line = process.stderr.readline()
            if line == '' and process.poll() != None:
                break
            if "worker" in line:
                break
            if "pestat" not in line:
                continue
            if "candidate unique pairs for" in line:
                if seen_candidate_line:
                    break
                seen_candidate_line = True
                nb_pairs = parse_list(line,4)
                for i in xrange(4):
                    stats[['FF', 'FR', 'RF', 'RR'][i]] = { 'nb_pairs' : nb_pairs[i] }
            if "orientation" in line:
                orientation = line.strip().split()[-1].replace('.','')
            if "mem_pestat] mean and std.dev:" in line:
                mean, stdev = parse_list(line,2)
                stats[orientation]['mean'] = mean
                stats[orientation]['stdev'] = stdev
                if orientation == 'RR':
                    # stats are complete
                    break
            sys.stdout.write(line)
            sys.stdout.flush()
        if process.poll() is None:
            process.terminate()
   
        results = sorted(stats.items(), key = lambda x: x[1]['nb_pairs'], reverse=True)
        most_likely = results[0]
        mean = most_likely[1]['mean']
        stdev = most_likely[1]['stdev']
        print "Orientation", most_likely[0], "mean", mean, "stdev", stdev



if __name__=="__main__":
    print __doc__
