#!/usr/bin/env python 
"""
This program aims to locate the path for data storage and experiment runs for each organisms in a workflow.

There are some predefined paths are loaded in this program statically. 
"""

import os 
from collections import defaultdict

def make_org_db(org_name_file, data_path, exp_path):
    """
    function to collect details of each organism

    @args org_name_file: text file containing organisms name
    @type org_name_file: str 
        
        ovis_aries      SRR645881       ensembl 78
        capra_hircus      SRR645885       ensembl 78

    @args data_path: data file storage path 
    @type data_path: str 
    @args exp_path: experiment related file path  
    @type exp_path: str 
    """

    org_fasta_file = dict( A_carolinensis = '%s/A_carolinensis/STARgenome/Anolis_carolinensis.AnoCar2.0.69.stable_genome.fa' % data_path,
    M_mulatta = "%s/M_mulatta/STARgenome/ensembl_release-69.fa" % data_path,
    O_cuniculus = "%s/O_cuniculus/STARgenome/ensembl_release-69.fa" % data_path,
    M_gallopavo = "%s/M_gallopavo/STARgenome/ensembl_release-69.fa" % data_path, 
    B_anthracis = '%s/B_anthracis/ensembl_release-21/Bacillus_anthracis_str_a0193.GCA_000181915.1.21.dna.toplevel.fa' % data_path,
    C_familiaris = '%s/C_familiaris/ensembl_release-69/Canis_familiaris.CanFam3.1.69.dna.toplevel.fa' % data_path,
    D_melanogaster = '%s/D_melanogaster/ensembl_release-69/Drosophila_melanogaster.BDGP5.69.dna.toplevel.fa' % data_path,
    E_caballus = '%s/E_caballus/ensembl_release-69/Equus_caballus.EquCab2.69_stable.fa.bz2' % data_path,
    M_domestica = '%s/M_domestica/ensembl_release-69/Monodelphis_domestica.BROADO5.69.dna.toplevel.fa' % data_path,
    O_sativa = '%s/O_sativa/phytozome_v9.0/Osativa_204.fa' % data_path,
    A_gambiae = '%s/A_gambiae/ensembl_release-20/Anopheles_gambiae.AgamP3.20.dna.fa' % data_path,
    B_rapa = '%s/B_rapa/phytozome_v9.0/Brapa_197_stable.fa' % data_path,
    C_japonica = '%s/C_japonica/STARgenome/Caenorhabditis_japonica.C_japonica-7.0.1.22.dna_sm.stable.fa' % data_path,
    G_gallus = '%s/G_gallus/ensembl_release-69/Gallus_gallus.WASHUC2.69_stable.fa.bz2' % data_path,
    M_musculus = '%s/M_musculus/ensembl_release-69/Mus_musculus.GRCm38.69_stable.fa.bz2' % data_path,  
    V_vinifera = '%s/V_vinifera/STARgenome/Vvinifera_145.fa' % data_path,
    A_mellifera = '%s/A_mellifera/ensembl_release-20/apiMel3_ucsc_chrom.fasta' % data_path,
    B_taurus = '%s/B_taurus/STARgenome/ensembl_release-69.fa' % data_path,
    C_rubella = '%s/C_rubella/phytozome_v9.0/Crubella_183.fa.gz' % data_path,
    D_rerio = '%s/D_rerio/ensembl_release-69/Danio_rerio.Zv9.69.dna.toplevel.fa.bz2' % data_path,
    G_max = '%s/G_max/phytozome_v9.0/Gmax_189_filter.fa' % data_path,
    M_truncatula = '%s/M_truncatula/STARgenome/Mtruncatula_198.fa' % data_path,
    P_pacificus = '%s/P_pacificus/STARgenome/Pristionchus_pacificus.P_pacificus-5.0.22.dna_sm.stable.fa' % data_path,
    S_scrofa = '%s/S_scrofa/ensembl_release-69/Sus_scrofa.Sscrofa10.2.69.dna_filter.fasta' % data_path,
    X_tropicalis = '%s/X_tropicalis/ensembl_release-69/Xenopus_tropicalis.JGI_4.2.69.dna_stable.fasta.bz2' % data_path,
    C_sativus = '%s/C_sativus/phytozome_v9.0/Csativus_122_filtered.fa' % data_path,
    D_simulans = '%s/D_simulans/ensembl_release-22/Drosophila_simulans.WUGSC1.22.dna_sm.dna.fa' % data_path,
    H_sapiens = '%s/H_sapiens/hg19_bowtie/hg19.fa' % data_path,
    O_anatinus = '%s/O_anatinus/ensembl_release-69/Ornithorhynchus_anatinus.OANA5.69-filtered_dna.fa' % data_path,
    P_troglodytes = '%s/P_troglodytes/STARgenome/Pan_troglodytes.CHIMP2.1.4.69_stable.fa' % data_path,
    S_tuberosum = '%s/S_tuberosum/phytozome_v9.0/Stuberosum_206.fa' % data_path,
    Z_mays = '%s/Z_mays/phytozome_v9.0/Zmays_181.fa' % data_path,
    A_thaliana = '%s/A_thaliana/arabidopsis_tair10/sequences/TAIR9_chr_all.fas' % data_path,
    C_elegans = '%s/C_elegans/ensembl_release-69/Caenorhabditis_elegans.WBcel215.69.dna.toplevel.fa' % data_path,
    D_discoideum = '%s/D_discoideum/STARgenome/Dictyostelium_discoideum.dictybase.01.21.dna.toplevel.fa' % data_path,
    D_yakuba = '%s/D_yakuba/STARgenome/Drosophila_yakuba.dyak_r1.3_FB2008_07.22.dna_sm_stable.fasta' % data_path,
    O_latipes = '%s/O_latipes/STARgenome/ensembl_release-74.fa' % data_path,
    R_norvegicus = '%s/R_norvegicus/ensembl_release-69/Rattus_norvegicus.RGSC3.4.69.dna.toplevel.fa.bz2' % data_path,
    C_briggsae = '%s/C_briggsae/STARgenome/Caenorhabditis_briggsae.CB4.22.dna_sm_stable.fasta' % data_path,
    C_brenneri = '%s/C_brenneri/STARgenome/Caenorhabditis_brenneri.C_brenneri-6.0.1b.22.dna_sm_stable.fa' % data_path,
    C_remanei = '%s/C_remanei/STARgenome/Caenorhabditis_remanei.C_remanei-15.0.1.22.dna_sm_stable.fa' % data_path,
    D_pseudoobscura = '%s/D_pseudoobscura/STARgenome/Drosophila_pseudoobscura.HGSC2.22.dna_sm_stable.fasta' % data_path,
    T_pseudonana = '%s/T_pseudonana/Thaps3/Thaps3_chromosomes_assembly_chromosomes_repeatmasked.fasta' % data_path,
    T_nigroviridis = '%s/T_nigroviridis/ensembl_release-69/Tetraodon_nigroviridis.TETRAODON8.69.dna.toplevel.fa' % data_path
    )

    star_index_file = dict( A_carolinensis = '%s/A_carolinensis/STARgenome/' % data_path,
    M_mulatta = "%s/M_mulatta/STARgenome/" % data_path,
    O_cuniculus = "%s/O_cuniculus/STARgenome/" % data_path,
    M_gallopavo = "%s/M_gallopavo/STARgenome/" % data_path, 
    B_anthracis = '%s/B_anthracis/STARgenome/' % data_path,
    C_familiaris = '%s/C_familiaris/STARgenome/' % data_path,
    D_melanogaster = '%s/D_melanogaster/STARgenome/' % data_path,
    E_caballus = '%s/E_caballus/STARgenome/' % data_path,
    M_domestica = '%s/M_domestica/STARgenome/' % data_path,
    O_sativa = '%s/O_sativa/STARgenome/' % data_path,
    A_gambiae = '%s/A_gambiae/STARgenome/' % data_path,
    B_rapa = '%s/B_rapa/STARgenome/' % data_path,
    C_japonica = '%s/C_japonica/STARgenome/' % data_path,
    G_gallus = '%s/G_gallus/STARgenome/' % data_path,
    M_musculus = '%s/M_musculus/ensembl_release-69/STARgenome' % data_path,  
    V_vinifera = '%s/V_vinifera/STARgenome/' % data_path,
    A_mellifera = '%s/A_mellifera/STARgenome/' % data_path,
    B_taurus = '%s/B_taurus/STARgenome/' % data_path,
    C_rubella = '%s/C_rubella/STARgenome/' % data_path,
    D_rerio = '%s/D_rerio/STARgenome/' % data_path,
    G_max = '%s/G_max/STARgenome/' % data_path,
    M_truncatula = '%s/M_truncatula/STARgenome/' % data_path,
    P_pacificus = '%s/P_pacificus/STARgenome/' % data_path,
    S_scrofa = '%s/S_scrofa/STARgenome/' % data_path,
    X_tropicalis = '%s/X_tropicalis/STARgenome/' % data_path,
    C_sativus = '%s/C_sativus/STARgenome/' % data_path,
    D_simulans = '%s/D_simulans/STARgenome/' % data_path,
    H_sapiens = '%s/H_sapiens/STARgenomes/hg19/' % data_path,
    O_anatinus = '%s/O_anatinus/STARgenome/' % data_path,
    P_troglodytes = '%s/P_troglodytes/STARgenome/' % data_path,
    S_tuberosum = '%s/S_tuberosum/STARgenome/' % data_path,
    Z_mays = '%s/Z_mays/STARgenome/' % data_path,
    A_thaliana = '%s/A_thaliana/STARgenome/' % data_path,
    C_elegans = '%s/C_elegans/STARgenome/' % data_path,
    D_discoideum = '%s/D_discoideum/STARgenome/' % data_path,
    D_yakuba = '%s/D_yakuba/STARgenome/' % data_path,
    O_latipes = '%s/O_latipes/STARgenome/' % data_path,
    R_norvegicus = '%s/R_norvegicus/STARgenome/' % data_path,
    C_briggsae = '%s/C_briggsae/STARgenome/' % data_path,
    C_brenneri = '%s/C_brenneri/STARgenome/' % data_path,
    C_remanei = '%s/C_remanei/STARgenome/' % data_path,
    D_pseudoobscura = '%s/D_pseudoobscura/STARgenome/' % data_path,
    T_pseudonana = '%s/T_pseudonana/STARgenome/' % data_path,
    T_nigroviridis = '%s/T_nigroviridis/STARgenome/' % data_path
    )

    org_gio_file = dict( A_carolinensis = '%s/A_carolinensis/ensembl_release-69/ens_gio/genome.config' % data_path,
    M_mulatta = "%s/M_mulatta/ensembl_release-69/ens_gio/genome.config" % data_path,
    O_cuniculus = "%s/O_cuniculus/ensembl_release-69/ens_gio/genome.config" % data_path,
    M_gallopavo = "%s/M_gallopavo/ensembl_release-69/ens_gio/genome.config" % data_path, 
    B_anthracis = '%s/B_anthracis/ensembl_release-21/ens_gio/genome.config' % data_path,
    C_familiaris = '%s/C_familiaris/ensembl_release-69/ens_gio/genome.config' % data_path,
    D_melanogaster = '%s/D_melanogaster/ensembl_release-69/ens_gio/genome.config' % data_path,
    E_caballus = '%s/E_caballus/ensembl_release-69/ens_gio/genome.config' % data_path,
    M_domestica = '%s/M_domestica/ensembl_release-69/ens_gio/genome.config' % data_path,
    O_sativa = '%s/O_sativa/phytozome_v9.0/phy_gio/genome.config' % data_path,
    A_gambiae = '%s/A_gambiae/ensembl_release-20/ens_gio/genome.config' % data_path,
    B_rapa = '%s/B_rapa/phytozome_v9.0/phy_gio/genome.config' % data_path,
    C_japonica = '%s/C_japonica/ensembl_release-22/ens_gio/genome.config' % data_path,
    G_gallus = '%s/G_gallus/ensembl_release-69/ens_gio/genome.config' % data_path,
    M_musculus = '%s/M_musculus/ensembl_release-69/ens_gio/genome.config' % data_path,  
    V_vinifera = '%s/V_vinifera/phytozome_v9.0/phy_gio/genome.config' % data_path,
    A_mellifera = '%s/A_mellifera/apiMel3_ucsc/ens_gio/genome.config' % data_path,
    B_taurus = '%s/B_taurus/ensembl_release-69/ens_gio/genome.config' % data_path,
    C_rubella = '%s/C_rubella/phytozome_v9.0/phy_gio/genome.config' % data_path,
    D_rerio = '%s/D_rerio/ensembl_release-69/ens_gio/genome.config' % data_path,
    G_max = '%s/G_max/phytozome_v9.0/phy_gio/genome.config' % data_path,
    M_truncatula = '%s/M_truncatula/ensembl_release-69/ens_gio/genome.config' % data_path,
    P_pacificus = '%s/P_pacificus/ensembl_release-22/ens_gio/genome.config' % data_path,
    S_scrofa = '%s/S_scrofa/ensembl_release-69/ens_gio/genome.config' % data_path,
    X_tropicalis = '%s/X_tropicalis/JGIv4-1/jgi_gio/genome.config' % data_path,
    C_sativus = '%s/C_sativus/phytozome_v9.0/phy_gio/genome.config' % data_path,
    D_simulans = '%s/D_simulans/ensembl_release-22/ens_gio/genome.config' % data_path,
    H_sapiens = '%s/H_sapiens/hg19_gio/genome.config' % data_path,
    O_anatinus = '%s/O_anatinus/ensembl_release-69/ens_gio/genome.config' % data_path,
    P_troglodytes = '%s/P_troglodytes/ensembl_release-69/ens_gio/genome.config' % data_path,
    S_tuberosum = '%s/S_tuberosum/phytozome_v9.0/S_tuberosum_gio/genome.config' % data_path,
    Z_mays = '%s/Z_mays/phytozome_v9.0/phy_gio/genome.config' % data_path,
    A_thaliana = '%s/A_thaliana/Tair10_gio/genome.config' % data_path,
    C_elegans = '%s/C_elegans/ensembl_release-69/ens_gio/genome.config' % data_path,
    D_discoideum = '%s/D_discoideum/ensembl_release-22/ens_gio/genome.config' % data_path,
    D_yakuba = '%s/D_yakuba/ensembl_release-22/ens_gio/genome.config' % data_path,
    O_latipes = '%s/O_latipes/ensembl_release-74/ens_gio/genome.config' % data_path,
    R_norvegicus = '%s/R_norvegicus/ensembl_release-69/ens_gio/genome.config' % data_path,
    C_briggsae = '%s/C_briggsae/ensembl_release-22/ens_gio/genome.config' % data_path,
    C_brenneri = '%s/C_brenneri/ensembl_release-22/ens_gio/genome.config' % data_path,
    C_remanei = '%s/C_remanei/ensembl_release-22/ens_gio/genome.config' % data_path,
    D_pseudoobscura = '%s/D_pseudoobscura/ensembl_release-22/ens_gio/genome.config' % data_path,
    T_pseudonana = '%s/T_pseudonana/Thaps3/Thaps3_gio/genome.config' % data_path,
    T_nigroviridis = '%s/T_nigroviridis/ensembl_release-69/ens_gio/genome.config' % data_path
    )

    org_gtf_file = dict( A_carolinensis = '%s/A_carolinensis/STARgenome/Anolis_carolinensis.AnoCar2.0.69.stable.gtf' % data_path,
    M_mulatta = "%s/M_mulatta/STARgenome/ensembl_release-69.gtf" % data_path,
    O_cuniculus = "%s/O_cuniculus/STARgenome/ensembl_release-69.gtf" % data_path,
    M_gallopavo = "%s/M_gallopavo/STARgenome/ensembl_release-69.gtf" % data_path, 
    B_anthracis = '%s/B_anthracis/ensembl_release-21/Bacillus_anthracis' % data_path,
    C_familiaris = '%s/C_familiaris/ensembl_release-69/Canis_familiaris.CanFam3.1.69.gtf' % data_path,
    D_melanogaster = '%s/D_melanogaster/ensembl_release-69/Drosophila_melanogaster.BDGP5.69.gtf' % data_path,
    E_caballus = '%s/E_caballus/ensembl_release-69/Equus_caballus.EquCab2.69.gtf' % data_path,
    M_domestica = '%s/M_domestica/ensembl_release-69/Monodelphis_domestica.BROADO5.69.gtf' % data_path,
    O_sativa = '%s/O_sativa/phytozome_v9.0/Osativa_204_gene.gff3' % data_path,
    A_gambiae = '%s/A_gambiae/ensembl_release-20/Anopheles_gambiae.AgamP3.20.gtf' % data_path,
    B_rapa = '%s/B_rapa/phytozome_v9.0/Brapa_197_gene.gff3' % data_path,
    C_japonica = '%s/C_japonica/ensembl_release-22/Caenorhabditis_japonica.C_japonica-7.0.1.22.gff3' % data_path,
    G_gallus = '%s/G_gallus/ensembl_release-69/Gallus_gallus.WASHUC2.69.gtf' % data_path,
    M_musculus = '%s/M_musculus/ensembl_release-69/Mus_musculus.GRCm38.69.gtf' % data_path,  
    V_vinifera = '%s/V_vinifera/STARgenome/Vvinifera_145_gene.gff3' % data_path,
    A_mellifera = '%s/A_mellifera/ensembl_release-20/apiMel3_ucsc_ensembl_genes.gtf' % data_path,
    B_taurus = '%s/ensembl_release-69/Bos_taurus.UMD3.1.69.gtf' % data_path,
    C_rubella = '%s/C_rubella/phytozome_v9.0/Crubella_183.gff3' % data_path,
    D_rerio = '%s/D_rerio/ensembl_release-69/Danio_rerio.Zv9.69.gtf' % data_path,
    G_max = '%s/G_max/phytozome_v9.0/Gmax_189_gene.gff3' % data_path,
    M_truncatula = '%s/M_truncatula/' % data_path,
    P_pacificus = '%s/P_pacificus/ensembl_release-22/Pristionchus_pacificus.P_pacificus-5.0.22.gtf' % data_path,
    S_scrofa = '%s/S_scrofa/ensembl_release-69/Sus_scrofa.Sscrofa10.2.69.gtf' % data_path,
    X_tropicalis = '%s/X_tropicalis/STARgenome/JGIv4-1.gtf' % data_path,
    C_sativus = '%s/C_sativus/phytozome_v9.0/Csativus_122_gene.gff3' % data_path,
    D_simulans = '%s/D_simulans/ensembl_release-22/Drosophila_simulans.WUGSC1.22.gff3' % data_path,
    H_sapiens = '%s/H_sapiens/ensembl_release-69/Homo_sapiens.GRCh37.69.gtf' % data_path,
    O_anatinus = '%s/O_anatinus/ensembl_release-69/' % data_path,
    P_troglodytes = '%s/P_troglodytes/ensembl_release-69/Pan_troglodytes.CHIMP2.1.4.69.gtf' % data_path,
    S_tuberosum = '%s/S_tuberosum/phytozome_v9.0/Stuberosum_206_gene.gff3' % data_path,
    Z_mays = '%s/Z_mays/phytozome_v9.0/Zmays_181_gene.gff3' % data_path,
    A_thaliana = '%s/A_thaliana/arabidopsis_tair10/annotations/TAIR10_GFF3_genes.gff' % data_path,
    C_elegans = '%s/C_elegans/ensembl_release-69/Caenorhabditis_elegans.WBcel215.69.gtf' % data_path,
    D_discoideum = '%s/D_discoideum/' % data_path,
    D_yakuba = '%s/D_yakuba/ensembl_release-22/Drosophila_yakuba.dyak_r1.3_FB2008_07.22.gff3' % data_path,
    O_latipes = '%s/O_latipes/STARgenome/Oryzias_latipes.MEDAKA1.74.gtf' % data_path,
    R_norvegicus = '%s/R_norvegicus/ensembl_release-69/Rattus_norvegicus.RGSC3.4.69.gtf' % data_path,
    C_briggsae = '%s/C_briggsae/ensembl_release-22/Caenorhabditis_briggsae.CB4.22.gff3' % data_path,
    C_brenneri = '%s/C_brenneri/ensembl_release-22/Caenorhabditis_brenneri.C_brenneri-6.0.1b.22.gff3' % data_path,
    C_remanei = '%s/C_remanei/ensembl_release-22/Caenorhabditis_remanei.C_remanei-15.0.1.22.gff3' % data_path,
    D_pseudoobscura = '%s/D_pseudoobscura/ensembl_release-22/Drosophila_pseudoobscura.HGSC2.22.gff3' % data_path,
    T_pseudonana = '%s/T_pseudonana/Thaps3/Thaps3_chromosomes_geneModels_FilteredModels2.gff' % data_path,
    T_nigroviridis = '%s/T_nigroviridis/ensembl_release-69/Tetraodon_nigroviridis.TETRAODON8.69.gtf' % data_path
    )

    org_db = defaultdict()

    ## get the organisms name and details on the experiment
    with open(org_name_file, "rU") as fh:
        for name in fh:
            name = name.strip('\n\r').split('\t')
            #print name 
            
            ## name shortening 
            token = name[0].split("_") 
            genus, species = token[0], token[-1]
            short_name = "%s_%s" % (genus[0].upper(), species.lower())  
             
            ## adding details 
            org_db[short_name] = dict(name = name[0])  
            org_db[short_name]['short_name'] = short_name

            ## sequencing reads files 
            sra_files = [] 
            if os.path.isdir("%s/%s/source_data" % (exp_path, short_name)):
                for sra_file in os.listdir("%s/%s/source_data" % (exp_path, short_name)):
                    file_prefix, ext = os.path.splitext(sra_file)
                    if ext == ".sra":
                        continue 
                    sra_files.append(sra_file) 
            else:
                # new organism, creating sub directories  
                for sub_dir in ['source_data', 'read_mapping', 'signal_labels', 'trans_pred']:
                    try:
                        os.makedirs("%s/%s/%s" % (exp_path, short_name, sub_dir))
                    except OSError:
                        print "Skipping creation of %s/%s/%s because it exists already." % (exp_path, short_name, sub_dir) 
                
            org_db[short_name]['fastq_path'] = "%s/%s/source_data" % (exp_path, short_name)
            org_db[short_name]['fastq'] = sra_files

            org_db[short_name]['star_wd'] = "%s/%s/read_mapping" % (exp_path, short_name)
            org_db[short_name]['trsk_wd'] = "%s/%s/trans_pred" % (exp_path, short_name)
            org_db[short_name]['labels_wd'] = "%s/%s/signal_labels" % (exp_path, short_name)

            org_db[short_name]['bam'] = "%s/%s/read_mapping/unique_map.bam" % (exp_path, short_name)
            org_db[short_name]['pred_gff'] = "%s/%s/trans_pred/ss_filter_predgenes.gff" % (exp_path, short_name)

            ## check for the genome sequence file 
            if short_name in org_fasta_file:
                org_db[short_name]['fasta'] = org_fasta_file[short_name]
            else:
                if not os.path.isdir("%s" % data_path):
                    os.makedirs("%s" % data_path) 
                else:
                    print "Skipping creation of %s because it exists already." % data_path

                org_db[short_name]['fasta'] = "%s" % data_path

            ## check for the genome index file 
            if short_name in star_index_file:
                org_db[short_name]['index'] = star_index_file[short_name]
            else:
                if not os.path.isdir("%s" % data_path):
                    os.makedirs("%s" % data_path)
                else:    
                    print "Skipping creation of %s because it exists already." % data_path

                org_db[short_name]['index'] = "%s/%s/STARgenome/" % (data_path, short_name) 

            ## TODO remove the dependency of gio from TSKM 
            if short_name in org_gio_file:
                org_db[short_name]['gio'] = org_gio_file[short_name]
            else:
                org_db[short_name]['gio'] = "%s" % data_path

            ## check the genome annotation 
            if short_name in org_gtf_file:
                org_db[short_name]['gtf'] = org_gtf_file[short_name]
                ## get the gtf feature lengths 
                if os.path.isfile(org_gtf_file[short_name]):
                    from fetch_remote_data import prepare_data as pd

                    feat_len_db = pd.make_anno_db(org_gtf_file[short_name]) 
                    org_db[short_name]['max_intron'] = feat_len_db['max_intron']
                    org_db[short_name]['max_exon'] = feat_len_db['max_exon']
            else:
                org_db[short_name]['gtf'] = data_path
                org_db[short_name]['max_intron'] = None
                org_db[short_name]['max_exon'] = None
            
            ## SRA/ENA run id 
            try:
                org_db[short_name]['sra_run_id'] = name[1]
            except:
                org_db[short_name]['sra_run_id'] = None 
                print "SRA run_id missing"

            ## genome annotation release number
            try:
                version = name[2].split(' ')
                org_db[short_name]['release_db'] = version[0]
                org_db[short_name]['release_num'] = version[-1]
                
                #sub_genome_folder = '%s/%s/%s_release_%s' % (data_path, short_name, version[0], version[-1]) 
                ## updating the fasta/gtf/index file location 
                #try:
                #    os.makedirs(sub_genome_folder)
                #    org_db[short_name]['fasta'] = sub_genome_folder
                #    org_db[short_name]['gtf'] = sub_genome_folder
                #    star_index_folder = '%s/STARgenome' % sub_genome_folder
                #    try:
                #        os.makedirs(star_index_folder)
                #        org_db[short_name]['index'] = star_index_folder
                #    except OSError:
                #        print "skipping creation of %s star index dir" % star_index_folder

                #except OSError:
                #    print "skipping creation of %s genome version dir" % sub_genome_folder

            except:
                org_db[short_name]['release_db'] = None
                org_db[short_name]['release_num'] = None 
                print "Genome annotation release database and number are missing"
                
    fh.close() 
    
    return org_db
