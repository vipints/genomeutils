#!/usr/bin/env python 
"""
This program aims to locate the path for data storage and experiment runs for each organisms in a workflow.
FIXME descriptions 
"""

import os 
import re 
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

    org_fasta_file = dict( A_carolinensis = '%s/A_carolinensis/ensembl_release-69/Anolis_carolinensis.AnoCar2.0.69.stable_genome.fa' % data_path,
    M_mulatta = "%s/M_mulatta/ensembl_release-69/ensembl_release-69.fa" % data_path,
    O_cuniculus = "%s/O_cuniculus/ensembl_release-69/ensembl_release-69.fa" % data_path,
    M_gallopavo = "%s/M_gallopavo/ensembl_release-69/ensembl_release-69.fa" % data_path, 
    B_anthracis = '%s/B_anthracis/ensembl_release-21/Bacillus_anthracis_str_a0193.GCA_000181915.1.21.dna.toplevel.fa' % data_path,
    C_familiaris = '%s/C_familiaris/ensembl_release-69/ensembl_release-69.fa' % data_path,
    D_melanogaster = '%s/D_melanogaster/ensembl_release-69/Drosophila_melanogaster.BDGP5.69.dna.toplevel.fa' % data_path,
    E_caballus = '%s/E_caballus/ensembl_release-69/Equus_caballus.EquCab2.69_stable.fa' % data_path,
    M_domestica = '%s/M_domestica/ensembl_release-69/Monodelphis_domestica.BROADO5.69.dna.toplevel.fa' % data_path,
    O_sativa = '%s/O_sativa/phytozome_v9.0/Osativa_204.fa' % data_path,
    A_gambiae = '%s/A_gambiae/anoGam1_ucsc/anoGam1_rm.fasta' % data_path,
    B_rapa = '%s/B_rapa/phytozome_v9.0/Brapa_197_stable.fa' % data_path,
    C_japonica = '%s/C_japonica/STARgenome/Caenorhabditis_japonica.C_japonica-7.0.1.22.dna_sm.stable.fa' % data_path,
    G_gallus = '%s/G_gallus/ensembl_release-69/Gallus_gallus.WASHUC2.69_stable.fa' % data_path,
    M_musculus = '%s/M_musculus/ensembl_release-69/Mus_musculus.GRCm38.69_stable.fa' % data_path,  
    V_vinifera = '%s/V_vinifera/phytozome_v9.0/Vvinifera_145.fa' % data_path,
    A_mellifera = '%s/A_mellifera/apiMel3_ucsc/apiMel3_sm.fasta' % data_path,
    B_taurus = '%s/B_taurus/ensembl_release-69/ensembl_release-69.fa' % data_path,
    C_rubella = '%s/C_rubella/phytozome_v9.0/Crubella_183.fa.gz' % data_path,
    D_rerio = '%s/D_rerio/ensembl_release-69/ensembl_release-69.fa' % data_path,
    G_max = '%s/G_max/phytozome_v9.0/Gmax_189_filter.fa' % data_path,
    M_truncatula = '%s/M_truncatula/STARgenome/Mtruncatula_198.fa' % data_path,
    P_pacificus = '%s/P_pacificus/STARgenome/Pristionchus_pacificus.P_pacificus-5.0.22.dna_sm.stable.fa' % data_path,
    S_scrofa = '%s/S_scrofa/ensembl_release-69/Sus_scrofa.Sscrofa10.2.69.dna_filter.fasta' % data_path,
    X_tropicalis = '%s/X_tropicalis/JGIv4-1/JGIv4-1.fa' % data_path,
    C_sativus = '%s/C_sativus/phytozome_v9.0/Csativus_122_filtered.fa' % data_path,
    D_simulans = '%s/D_simulans/ensembl_release-22/Drosophila_simulans.WUGSC1.22.dna_sm.dna.fa' % data_path,
    H_sapiens = '%s/H_sapiens/hg19_bowtie/hg19.fa' % data_path,
    O_anatinus = '%s/O_anatinus/ensembl_release-69/Ornithorhynchus_anatinus.OANA5.69-filtered_dna.fa' % data_path,
    N_vitripennis = '%s/N_vitripennis/ensembl_release_22/N_vitripennis_dna_sm.fa' % data_path,
    P_troglodytes = '%s/P_troglodytes/ensembl_release-69/Pan_troglodytes.CHIMP2.1.4.69_stable.fa' % data_path,
    S_tuberosum = '%s/S_tuberosum/phytozome_v9.0/Stuberosum_206.fa' % data_path,
    Z_mays = '%s/Z_mays/phytozome_v9.0/Zmays_181.fa' % data_path,
    A_thaliana = '%s/A_thaliana/arabidopsis_tair10/sequences/TAIR9_chr_all.fas' % data_path,
    O_aries = '%s/O_aries/ensembl_release_78/O_aries_dna_sm.fa' % data_path,
    C_jacchus = '%s/C_jacchus/ensembl_release_69/C_jacchus_dna_sm.fa' % data_path,
    C_elegans = '%s/C_elegans/ensembl_release-69/Caenorhabditis_elegans.WBcel215.69.dna.toplevel.fa' % data_path,
    O_latipes = '%s/O_latipes/ensembl_release-74/ensembl_release-74.fa' % data_path,
    R_norvegicus = '%s/R_norvegicus/ensembl_release-69/Rattus_norvegicus.RGSC3.4.69.dna.toplevel.fa' % data_path,
    C_briggsae = '%s/C_briggsae/ensembl_release-22/Caenorhabditis_briggsae.CB4.22.dna_sm_stable.fasta' % data_path,
    T_nigroviridis = '%s/T_nigroviridis/ensembl_release-69/Tetraodon_nigroviridis.TETRAODON8.69.dna.toplevel.fa' % data_path
    )

    org_gtf_file = dict( A_carolinensis = '%s/A_carolinensis/ensembl_release-69/Anolis_carolinensis.AnoCar2.0.69.stable.gtf' % data_path,
    M_mulatta = "%s/M_mulatta/ensembl_release-69/ensembl_release-69.gtf" % data_path,
    O_cuniculus = "%s/O_cuniculus/ensembl_release-69/ensembl_release-69.gtf" % data_path,
    M_gallopavo = "%s/M_gallopavo/ensembl_release-69/ensembl_release-69.gtf" % data_path, 
    B_anthracis = '%s/B_anthracis/ensembl_release-21/Bacillus_anthracis' % data_path,
    C_familiaris = '%s/C_familiaris/ensembl_release-69/Canis_familiaris.CanFam3.1.69.gtf' % data_path,
    D_melanogaster = '%s/D_melanogaster/ensembl_release-69/Drosophila_melanogaster.BDGP5.69.gtf' % data_path,
    E_caballus = '%s/E_caballus/ensembl_release-69/Equus_caballus.EquCab2.69.gtf' % data_path,
    M_domestica = '%s/M_domestica/ensembl_release-69/Monodelphis_domestica.BROADO5.69.gtf' % data_path,
    O_sativa = '%s/O_sativa/phytozome_v9.0/Osativa_204_gene.gff3' % data_path,
    A_gambiae = '%s/A_gambiae/anoGam1_ucsc/anoGam1_ucsc.gtf' % data_path,
    B_rapa = '%s/B_rapa/phytozome_v9.0/Brapa_197_gene.gff3' % data_path,
    C_japonica = '%s/C_japonica/ensembl_release-22/Caenorhabditis_japonica.C_japonica-7.0.1.22.gff3' % data_path,
    G_gallus = '%s/G_gallus/ensembl_release-69/Gallus_gallus.WASHUC2.69.gtf' % data_path,
    M_musculus = '%s/M_musculus/ensembl_release-69/Mus_musculus.GRCm38.69.gtf' % data_path,  
    V_vinifera = '%s/V_vinifera/phytozome_v9.0/Vvinifera_145_gene.gff3' % data_path,
    A_mellifera = '%s/A_mellifera/apiMel3_ucsc/apiMel2_ucsc.gtf' % data_path,
    B_taurus = '%s/B_taurus/ensembl_release-69/Bos_taurus.UMD3.1.69.gtf' % data_path,
    C_jacchus = '%s/C_jacchus/ensembl_release_69/C_jacchus.gff' % data_path,
    C_rubella = '%s/C_rubella/phytozome_v9.0/Crubella_183.gff3' % data_path,
    D_rerio = '%s/D_rerio/ensembl_release-69/ensembl_release-69.gtf' % data_path,
    G_max = '%s/G_max/phytozome_v9.0/Gmax_189_gene.gff3' % data_path,
    N_vitripennis = '%s/N_vitripennis/ensembl_release_22/N_vitripennis.gtf' % data_path,
    O_aries = '%s/O_aries/ensembl_release_78/O_aries.gtf' % data_path,
    M_truncatula = '%s/M_truncatula/' % data_path,
    P_pacificus = '%s/P_pacificus/ensembl_release-22/Pristionchus_pacificus.P_pacificus-5.0.22.gtf' % data_path,
    S_scrofa = '%s/S_scrofa/ensembl_release-69/Sus_scrofa.Sscrofa10.2.69.gtf' % data_path,
    X_tropicalis = '%s/X_tropicalis/JGIv4-1/JGIv4-1.gff' % data_path,
    C_sativus = '%s/C_sativus/phytozome_v9.0/Csativus_122_gene.gff3' % data_path,
    D_simulans = '%s/D_simulans/ensembl_release-22/Drosophila_simulans.WUGSC1.22.gff3' % data_path,
    H_sapiens = '%s/H_sapiens/ensembl_release-69/Homo_sapiens.GRCh37.69_stable.gtf' % data_path,
    O_anatinus = '%s/O_anatinus/ensembl_release-69/' % data_path,
    P_troglodytes = '%s/P_troglodytes/ensembl_release-69/Pan_troglodytes.CHIMP2.1.4.69.gtf' % data_path,
    S_tuberosum = '%s/S_tuberosum/phytozome_v9.0/Stuberosum_206_gene.gff3' % data_path,
    Z_mays = '%s/Z_mays/phytozome_v9.0/Zmays_181_gene.gff3' % data_path,
    A_thaliana = '%s/A_thaliana/arabidopsis_tair10/annotations/TAIR10_GFF3_genes.gff' % data_path,
    C_elegans = '%s/C_elegans/ensembl_release-69/Caenorhabditis_elegans.WBcel215.69.gtf' % data_path,
    D_discoideum = '%s/D_discoideum/' % data_path,
    D_yakuba = '%s/D_yakuba/ensembl_release-22/Drosophila_yakuba.dyak_r1.3_FB2008_07.22.gff3' % data_path,
    O_latipes = '%s/O_latipes/ensembl_release-74/Oryzias_latipes.MEDAKA1.74.gtf' % data_path,
    R_norvegicus = '%s/R_norvegicus/ensembl_release-69/Rattus_norvegicus.RGSC3.4.69.gtf' % data_path,
    C_briggsae = '%s/C_briggsae/ensembl_release-22/Caenorhabditis_briggsae.CB4.22.gff3' % data_path,
    C_brenneri = '%s/C_brenneri/ensembl_release-22/Caenorhabditis_brenneri.C_brenneri-6.0.1b.22.gff3' % data_path,
    C_remanei = '%s/C_remanei/ensembl_release-22/Caenorhabditis_remanei.C_remanei-15.0.1.22.gff3' % data_path,
    D_pseudoobscura = '%s/D_pseudoobscura/ensembl_release-22/Drosophila_pseudoobscura.HGSC2.22.gff3' % data_path,
    T_pseudonana = '%s/T_pseudonana/Thaps3/Thaps3_chromosomes_geneModels_FilteredModels2.gff' % data_path,
    T_nigroviridis = '%s/T_nigroviridis/ensembl_release-69/Tetraodon_nigroviridis.TETRAODON8.69.gtf' % data_path
    )

    ## TODO algorithms details 

    ## experiment details  
    org_db = defaultdict()
    
    for ent in config_map['experiment']:
        short_name = ent['organism_name'] 
        sra_run_id = ent['sra_run_id']
        genome_build_version = ent['genome_build_version']

        org_db[short_name] = dict(short_name = short_name)  
        org_db[short_name]['sra_run_id'] = sra_run_id
        org_db[short_name]['genome_release_db'] = genome_build_version

        ## sequencing reads files 
        sra_files = [] 
        if os.path.isdir("%s/%s/source_data" % (exp_path, short_name)):
            for sra_file in os.listdir("%s/%s/source_data" % (exp_path, short_name)):
                file_prefx, ext = os.path.splitext(sra_file)
                if ext == ".sra": ## skipping the original .sra binary file 
                    continue
                if re.search(sra_run_id, sra_file):
                    sra_files.append(sra_file) 
        else:
            print "warning: didn't find sequencing read files %s/%s/source_data" % (exp_path, short_name) 
                
        org_db[short_name]['fastq_path'] = "%s/%s/source_data" % (exp_path, short_name)
        org_db[short_name]['fastq'] = sra_files

        ## read mapping, read assembly and label generation working folders 
        org_db[short_name]['read_map_dir'] = "%s/%s/read_mapping" % (exp_path, short_name)
        org_db[short_name]['read_assembly_dir'] = "%s/%s/trans_pred" % (exp_path, short_name)
        org_db[short_name]['labels_dir'] = "%s/%s/signal_labels" % (exp_path, short_name)

        ## calculate the sequence read length
        readlength = 0 
        if opt_action in ["c", "a", "2", "3"]: ## perform this action only for selected options 
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
            org_db[short_name]['fasta'] = org_fasta_file[short_name]
        else:
            print "warning: didn't find genome sequence file fasta/fa for organism %s" % short_name
            org_db[short_name]['fasta'] = "%s/%s/%s" % (data_path, short_name, genome_build_version) 

        if not os.path.isdir("%s/%s/%s/STARgenome" % (data_path, short_name, genome_build_version)):
            os.makedirs("%s/%s/%s/STARgenome" % (data_path, short_name, genome_build_version))

        org_db[short_name]['genome_index_dir'] = "%s/%s/%s/STARgenome/" % (data_path, short_name, genome_build_version) 

        ## check the genome annotation 
        if short_name in org_gtf_file:
            org_db[short_name]['gtf'] = org_gtf_file[short_name]

            if opt_action in ["c", "a", "2", "3"]: ## perform this action only for selected options 
                ## get the gtf feature lengths 
                if os.path.isfile(org_gtf_file[short_name]):
                    from fetch_remote_data import prepare_data as pd
                    feat_len_db = pd.make_anno_db(org_gtf_file[short_name]) 
                    org_db[short_name]['max_intron_len'] = feat_len_db['max_intron']
                    org_db[short_name]['max_exon_len'] = feat_len_db['max_exon']
                else:
                    print "error: the provided gtf file %s is not available to read. Please check!" % org_gtf_file[short_name]
                    sys.exit(-1)
        else:
                print "warning: didn't find an annotation file gtf/gff for organism %s" % short_name
                org_db[short_name]['gtf'] = "%s/%s/%s" % (data_path, short_name, genome_build_version)
                org_db[short_name]['max_intron_len'] = None
                org_db[short_name]['max_exon_len'] = None
        
        print "fetched details for %s" % short_name 
            
    return org_db
