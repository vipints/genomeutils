#!/usr/bin/env python 
"""
This program contains some predefined variables which hold the relative path for 
data storage and experiment runs for each organisms in a workflow. The idea behind 
this is to have a speed query to the organisms detailed informations, specifically 
genome annotation file etc.

The program includes 
    - counting the sequencing read length for the provided experiment
    - getting the details of annotated features 
    - preparing for the experiment relative paths

Requires:
    biopython 
    gfftools 
"""

import os 
import re 
import sys 
import yaml 
from Bio import SeqIO
from gfftools import helper 
from collections import defaultdict

def experiment_db(config_file, opt_action):
    """
    function to collect details of each organism

    FIXME descriptions 
    @args config_file: yaml file contain the information for the experiment
    @type config_file: str 
    """

    ## parsing the config file 
    config_map = yaml.safe_load(open(config_file, "rU"))

    data_path = config_map['genome_data_path']['dir']
    exp_path = config_map['experiment_data_path']['dir']

    org_fasta_file = dict( A_carolinensis = '%s/A_carolinensis/ensembl_release_79/ensembl_release_79.fas' % data_path,
    M_mulatta = "%s/M_mulatta/ensembl_release_79/ensembl_release_79.fas" % data_path,
    O_cuniculus = "%s/O_cuniculus/ensembl_release_79/ensembl_release_79.fas" % data_path,
    M_gallopavo = "%s/M_gallopavo/ensembl_release_79/ensembl_release_79.fas.bz2" % data_path, 
    B_anthracis = '%s/B_anthracis/ensembl_release-21/Bacillus_anthracis_str_a0193.GCA_000181915.1.21.dna.toplevel.fa' % data_path,
    C_familiaris = '%s/C_familiaris/ensembl_release_79/ensembl_release_79.fas.bz2' % data_path,
    D_melanogaster = '%s/D_melanogaster/ensembl_release_79/ensembl_release_79.fas' % data_path,
    E_caballus = '%s/E_caballus/ensembl_release_79/ensembl_release_79.fas.bz2' % data_path,
    M_domestica = '%s/M_domestica/ensembl_release_79/ensembl_release_79.fas' % data_path,
    O_sativa = '%s/O_sativa/phytozome_v9.0/Osativa_204.fa' % data_path,
    A_gambiae = '%s/A_gambiae/ensembl_release_28/ensembl_release_28.fas' % data_path,
    B_rapa = '%s/B_rapa/phytozome_v9.0/Brapa_197_stable.fa' % data_path,
    G_gallus = '%s/G_gallus/ensembl_release_79/ensembl_release_79.fas' % data_path,
    M_musculus = '%s/M_musculus/ensembl_release_79/ensembl_release_79.fas' % data_path,  
    V_vinifera = '%s/V_vinifera/phytozome_v9.0/phytozome_v9.0.fas.bz2' % data_path,
    A_mellifera = '%s/A_mellifera/ensembl_release_28/ensembl_release_28.fas' % data_path,
    B_taurus = '%s/B_taurus/ensembl_release_79/ensembl_release_79.fas.bz2' % data_path,
    C_rubella = '%s/C_rubella/phytozome_v9.0/Crubella_183.fa.gz' % data_path,
    D_rerio = '%s/D_rerio/ensembl_release_79/ensembl_release_79.fas' % data_path,
    G_max = '%s/G_max/phytozome_v9.0/Gmax_189_filter.fa' % data_path,
    M_truncatula = '%s/M_truncatula/STARgenome/Mtruncatula_198.fa' % data_path,
    P_pacificus = '%s/P_pacificus/STARgenome/Pristionchus_pacificus.P_pacificus-5.0.22.dna_sm.stable.fa' % data_path,
    S_scrofa = '%s/S_scrofa/ensembl_release_79/ensembl_release_79.fas' % data_path,
    X_tropicalis = '%s/X_tropicalis/JGIv4-1/JGIv4-1.fa' % data_path,
    C_sativus = '%s/C_sativus/phytozome_v9.0/Csativus_122_filtered.fa' % data_path,
    D_simulans = '%s/D_simulans/ensembl_release_28/ensembl_release_28.fas' % data_path,
    H_sapiens = '%s/H_sapiens/hg19_bowtie2/hg19.fa' % data_path,
    N_vitripennis = '%s/N_vitripennis/ensembl_release_22/N_vitripennis_dna_sm.fa' % data_path,
    P_troglodytes = '%s/P_troglodytes/ensembl_release_79/ensembl_release_79.fas' % data_path,
    S_tuberosum = '%s/S_tuberosum/phytozome_v9.0/Stuberosum_206.fa' % data_path,
    Z_mays = '%s/Z_mays/phytozome_v9.0/Zmays_181.fa' % data_path,
    A_thaliana = '%s/A_thaliana/arabidopsis_tair10/sequences/TAIR9_chr_all.fas' % data_path,
    O_aries = '%s/O_aries/ensembl_release_79/ensembl_release_79.fas' % data_path,
    C_jacchus = '%s/C_jacchus/ensembl_release_79/ensembl_release_79.fas' % data_path,
    C_elegans = '%s/C_elegans/ensembl_release_79/ensembl_release_79.fas' % data_path,
    O_latipes = '%s/O_latipes/ensembl_release_79/ensembl_release_79.fas' % data_path,
    R_norvegicus = '%s/R_norvegicus/ensembl_release_79/ensembl_release_79.fas' % data_path,
    G_gorilla = '%s/G_gorilla/ensembl_release_79/ensembl_release_79.fas' % data_path,
    P_paniscus = '%s/P_paniscus/eva_mpg_de/eva_mpg_de.fas' % data_path,
    C_porcellus = '%s/C_porcellus/ensembl_release_79/ensembl_release_79.fas' % data_path,
    O_anatinus = '%s/O_anatinus/ensembl_release_79/ensembl_release_79.fas' % data_path,
    A_platyrhynchos = '%s/A_platyrhynchos/ensembl_release_79/ensembl_release_79.fas' % data_path,
    O_niloticus = '%s/O_niloticus/ensembl_release_79/ensembl_release_79.fas' % data_path,
    L_chalumnae = '%s/L_chalumnae/ensembl_release_79/ensembl_release_79.fas' % data_path,
    H_glaber = '%s/H_glaber/naked_mole_rat_db/naked_mole_rat_db.fas' % data_path,
    M_eugenii = '%s/M_eugenii/ensembl_release_79/ensembl_release_79.fas' % data_path,
    C_briggsae = '%s/C_briggsae/ensembl_release_28/ensembl_release_28.fas' % data_path,
    C_japonica = '%s/C_japonica/ensembl_release_28/ensembl_release_28.fas' % data_path,
    C_remanei = '%s/C_remanei/ensembl_release_28/ensembl_release_28.fas' % data_path,
    P_marinus = '%s/P_marinus/ensembl_release_79/ensembl_release_79.fas' % data_path,
    C_brenneri = '%s/C_brenneri/ensembl_release_28/ensembl_release_28.fas' % data_path,
    C_intestinalis = '%s/C_intestinalis/ensembl_release_79/ensembl_release_79.fas' % data_path,
    S_cerevisiae = '%s/S_cerevisiae/ensembl_release_79/ensembl_release_79.fas' % data_path,
    S_pombe = '%s/S_pombe/ensembl_release_28/ensembl_release_28.fas' % data_path,
    A_aegypti = '%s/A_aegypti/ensembl_release_28/ensembl_release_28.fas' % data_path,
    T_nigroviridis = '%s/T_nigroviridis/ensembl_release_79/ensembl_release_79.fas' % data_path
    )

    org_gtf_file = dict( A_carolinensis = '%s/A_carolinensis/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    M_mulatta = "%s/M_mulatta/ensembl_release_79/ensembl_release_79.gtf" % data_path,
    O_cuniculus = "%s/O_cuniculus/ensembl_release_79/ensembl_release_79.gtf" % data_path,
    M_gallopavo = "%s/M_gallopavo/ensembl_release_79/ensembl_release_79.gtf.bz2" % data_path, 
    B_anthracis = '%s/B_anthracis/ensembl_release-21/Bacillus_anthracis' % data_path,
    C_familiaris = '%s/C_familiaris/ensembl_release_79/ensembl_release_79.gtf.bz2' % data_path,
    D_melanogaster = '%s/D_melanogaster/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    E_caballus = '%s/E_caballus/ensembl_release_79/ensembl_release_79.gtf.bz2' % data_path,
    M_domestica = '%s/M_domestica/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    O_sativa = '%s/O_sativa/phytozome_v9.0/Osativa_204_gene.gff3' % data_path,
    A_gambiae = '%s/A_gambiae/ensembl_release_28/ensembl_release_28.gtf' % data_path,
    B_rapa = '%s/B_rapa/phytozome_v9.0/Brapa_197_gene.gff3' % data_path,
    G_gallus = '%s/G_gallus/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    M_musculus = '%s/M_musculus/ensembl_release_79/ensembl_release_79.gtf' % data_path,  
    V_vinifera = '%s/V_vinifera/phytozome_v9.0/phytozome_v9.0.gff.bz2' % data_path,
    A_mellifera = '%s/A_mellifera/ensembl_release_28/ensembl_release_28.gtf' % data_path,
    B_taurus = '%s/B_taurus/ensembl_release_79/ensembl_release_79.gtf.bz2' % data_path,
    C_jacchus = '%s/C_jacchus/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    C_rubella = '%s/C_rubella/phytozome_v9.0/Crubella_183.gff3' % data_path,
    D_rerio = '%s/D_rerio/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    G_max = '%s/G_max/phytozome_v9.0/Gmax_189_gene.gff3' % data_path,
    N_vitripennis = '%s/N_vitripennis/ensembl_release_22/N_vitripennis.gtf' % data_path,
    O_aries = '%s/O_aries/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    M_truncatula = '%s/M_truncatula/' % data_path,
    P_pacificus = '%s/P_pacificus/ensembl_release-22/Pristionchus_pacificus.P_pacificus-5.0.22.gtf' % data_path,
    S_scrofa = '%s/S_scrofa/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    X_tropicalis = '%s/X_tropicalis/JGIv4-1/JGIv4-1.gff' % data_path,
    C_sativus = '%s/C_sativus/phytozome_v9.0/Csativus_122_gene.gff3' % data_path,
    D_simulans = '%s/D_simulans/ensembl_release_28/ensembl_release_28.gtf' % data_path,
    H_sapiens = '%s/H_sapiens/ensembl_release_79/ensembl_release_79.gtf.bz2' % data_path,
    P_troglodytes = '%s/P_troglodytes/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    S_tuberosum = '%s/S_tuberosum/phytozome_v9.0/Stuberosum_206_gene.gff3' % data_path,
    Z_mays = '%s/Z_mays/phytozome_v9.0/Zmays_181_gene.gff3' % data_path,
    A_thaliana = '%s/A_thaliana/arabidopsis_tair10/annotations/TAIR10_GFF3_genes.gff' % data_path,
    C_elegans = '%s/C_elegans/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    D_discoideum = '%s/D_discoideum/' % data_path,
    D_yakuba = '%s/D_yakuba/ensembl_release-22/Drosophila_yakuba.dyak_r1.3_FB2008_07.22.gff3' % data_path,
    O_latipes = '%s/O_latipes/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    R_norvegicus = '%s/R_norvegicus/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    D_pseudoobscura = '%s/D_pseudoobscura/ensembl_release-22/Drosophila_pseudoobscura.HGSC2.22.gff3' % data_path,
    T_pseudonana = '%s/T_pseudonana/Thaps3/Thaps3_chromosomes_geneModels_FilteredModels2.gff' % data_path,
    G_gorilla = '%s/G_gorilla/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    C_porcellus = '%s/C_porcellus/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    O_anatinus = '%s/O_anatinus/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    A_platyrhynchos = '%s/A_platyrhynchos/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    O_niloticus = '%s/O_niloticus/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    L_chalumnae = '%s/L_chalumnae/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    C_briggsae = '%s/C_briggsae/ensembl_release_28/ensembl_release_28.gtf' % data_path,
    M_eugenii = '%s/M_eugenii/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    C_remanei = '%s/C_remanei/ensembl_release_28/ensembl_release_28.gtf' % data_path,
    C_brenneri = '%s/C_brenneri/ensembl_release_28/ensembl_release_28.gtf' % data_path,
    C_intestinalis = '%s/C_intestinalis/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    S_pombe = '%s/S_pombe/ensembl_release_28/ensembl_release_28.gtf' % data_path,
    P_marinus = '%s/P_marinus/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    S_cerevisiae = '%s/S_cerevisiae/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    C_japonica = '%s/C_japonica/ensembl_release_28/ensembl_release_28.gtf' % data_path,
    A_aegypti = '%s/A_aegypti/ensembl_release_28/ensembl_release_28.gtf' % data_path,
    T_nigroviridis = '%s/T_nigroviridis/ensembl_release_79/ensembl_release_79.gtf' % data_path
    )

    ## TODO algorithms details 
    #A_aegypti = '%s/A_aegypti/ensembl_release_28/ensembl_release_28.gtf' % data_path,
    #S_cerevisiae = '%s/S_cerevisiae/ensembl_release_79/ensembl_release_79.gtf' % data_path,
    #H_glaber = '%s/H_glaber/naked_mole_rat_db/naked_mole_rat_db.gtf' % data_path,
    #H_glaber = '%s/H_glaber/naked_mole_rat_db/naked_mole_rat_db.gtf' % data_path,

    ## experiment details  
    org_db = defaultdict()
    
    for ent in config_map['experiment']:
        species_name = ent['organism_name'] 
        sra_run_id = ent['sra_run_id']
        genome_build_version = ent['genome_build_version']
        db_server = ent['release_db']

        ## mapping to short names       arabidopsis_thaliana --> A_thaliana
        genus, species = species_name.strip().split("_")
        short_name = "%s_%s" % (genus[0].upper(), species)

        org_db[short_name] = dict(short_name = short_name)  
        org_db[short_name]['long_name'] = species_name
        org_db[short_name]['sra_run_id'] = sra_run_id
        org_db[short_name]['genome_release_db'] = genome_build_version
        ## the broad path to the experiment 
        org_db[short_name]['genome_dir'] = data_path
        org_db[short_name]['experiment_dir'] = exp_path

        build_release = genome_build_version.split("_")
        org_db[short_name]['release_db'] = db_server ## ensembl_metazoa, phytozome
        org_db[short_name]['release_nb'] = build_release[-1] ## build number 

        sra_files = [] ## sequencing reads files 
        if os.path.isdir("%s/%s/source_data" % (exp_path, short_name)):
            for sra_file in os.listdir("%s/%s/source_data" % (exp_path, short_name)):
                file_prefx, ext = os.path.splitext(sra_file)
                if ext == ".sra": ## skipping the original .sra binary file 
                    continue
                if re.search(sra_run_id, sra_file):
                    sra_files.append(sra_file) 
        else:
            print "warning: missing sequencing read trace file %s/%s/source_data" % (exp_path, short_name) 
                
        org_db[short_name]['fastq_path'] = "%s/%s/source_data" % (exp_path, short_name)
        org_db[short_name]['fastq'] = sra_files

        ## read mapping, read assembly and label generation working folders 
        for sub_dir in ['read_mapping', 'signal_labels', 'trans_pred', 'source_data']:
            work_path = "%s/%s/%s" % (exp_path, short_name, sub_dir)

            if not os.path.isdir(work_path):
                try:
                    os.makedirs(work_path)
                except OSError:
                    print "error: cannot create the directory %s." % work_path
                    sys.exit(0)

        org_db[short_name]['read_map_dir'] = "%s/%s/read_mapping" % (exp_path, short_name)
        org_db[short_name]['read_assembly_dir'] = "%s/%s/trans_pred" % (exp_path, short_name)
        org_db[short_name]['labels_dir'] = "%s/%s/signal_labels" % (exp_path, short_name)

        ## calculate the sequence read length
        readlength = 0 
        if opt_action in ["2", "3"]: ## perform this action only for selected options 
            if sra_files:
                fqfile = os.path.join(org_db[short_name]['fastq_path'], sra_files[0])
                print 'using sequencing read file %s to determine readLength' % fqfile
                fh = helper.open_file(fqfile)
                for rec in SeqIO.parse(fh, "fastq"):
                    readlength = len(rec.seq)
                    break
                fh.close() 
        org_db[short_name]['read_length'] = readlength

        ## check for the genome sequence file 
        if short_name in org_fasta_file:
            if os.path.isfile(org_fasta_file[short_name]):
                org_db[short_name]['fasta'] = org_fasta_file[short_name]
            else:
                org_db[short_name]['fasta'] = None
        else:
            print "warning: missing genome sequence file for %s under %s/%s/%s" % (short_name, data_path, short_name, genome_build_version)
            org_db[short_name]['fasta'] = None

        if not os.path.isdir("%s/%s/%s/STARgenome" % (data_path, short_name, genome_build_version)):
            try:
                os.makedirs("%s/%s/%s/STARgenome" % (data_path, short_name, genome_build_version))
            except OSError:
                print "error: cannot create the directory %s/%s/%s" % (data_path, short_name, genome_build_version) 
                sys.exit(0)

        org_db[short_name]['genome_index_dir'] = "%s/%s/%s/STARgenome/" % (data_path, short_name, genome_build_version) 

        ## check the genome annotation 
        if short_name in org_gtf_file:
            if os.path.isfile(org_gtf_file[short_name]):
                org_db[short_name]['gtf'] = org_gtf_file[short_name]

                if opt_action in ["2", "3", "4", "c"]: ## perform this action only for selected options 
                    ## get the gtf feature lengths 
                    from fetch_remote_data import prepare_data as pd
                    feat_len_db = pd.make_anno_db(org_gtf_file[short_name]) 
                    org_db[short_name]['max_intron_len'] = feat_len_db['max_intron']
                    org_db[short_name]['max_exon_len'] = feat_len_db['max_exon']
            else:
                exit("error: the provided gtf file %s is not available to read. Please check!" % org_gtf_file[short_name])
        else:
            print("warning: missing annotation file for %s under %s/%s/%s" % (short_name, data_path, short_name, genome_build_version))
            org_db[short_name]['gtf'] = None
            org_db[short_name]['max_intron_len'] = None
            org_db[short_name]['max_exon_len'] = None
        
        print("fetched details for %s" % short_name) 
            
    return org_db
