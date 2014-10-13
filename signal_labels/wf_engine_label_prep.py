#!/usr/bin/env python 
"""
Workflow script for generating different genomic signal labels.

Usage: python wf_engine_label_prep.py A_thaliana
"""

import sys 
from signal_labels import generate_genome_seq_labels

try:
    org_name = sys.argv[1]
except:
    
    print __doc__ 
    sys.exit(-1)

org_fasta_file = dict( A_carolinensis = '/cbio/grlab/share/databases/genomes/A_carolinensis/STARgenome/Anolis_carolinensis.AnoCar2.0.69.stable_genome.fa',
    M_mulatta = "/cbio/grlab/share/databases/genomes/M_mulatta/STARgenome/ensembl_release-69.fa",
    O_cuniculus = "/cbio/grlab/share/databases/genomes/O_cuniculus/STARgenome/ensembl_release-69.fa",
    M_gallopavo = "/cbio/grlab/share/databases/genomes/M_gallopavo/STARgenome/ensembl_release-69.fa",
    M_gallopavo = "/cbio/grlab/share/databases/genomes/M_gallopavo/STARgenome/ensembl_release-69.fa",
    B_anthracis = '/cbio/grlab/share/databases/genomes/B_anthracis/ensembl_release-21/Bacillus_anthracis_str_a0193.GCA_000181915.1.21.dna.toplevel.fa',
    C_familiaris = '/cbio/grlab/share/databases/genomes/C_familiaris/STARgenome/ensembl_release-69.fa',
    D_melanogaster = '/cbio/grlab/share/databases/genomes/D_melanogaster/ensembl_release-69/Drosophila_melanogaster.BDGP5.69.dna.toplevel.fa',
    E_caballus = '/cbio/grlab/share/databases/genomes/E_caballus/ensembl_release-69/Equus_caballus.EquCab2.69_stable.fa.bz2',
    M_domestica = '/cbio/grlab/share/databases/genomes/M_domestica/ensembl_release-69/Monodelphis_domestica.BROADO5.69.dna.toplevel.fa',
    O_sativa = '/cbio/grlab/share/databases/genomes/O_sativa/phytozome_v9.0/Osativa_204.fa',
    A_gambiae = '/cbio/grlab/share/databases/genomes/A_gambiae/ensembl_release-20/Anopheles_gambiae.AgamP3.20.dna.fa',
    B_rapa = '/cbio/grlab/share/databases/genomes/B_rapa/phytozome_v9.0/Brapa_197_stable.fa',
    C_japonica = '/cbio/grlab/share/databases/genomes/C_japonica/STARgenome/Caenorhabditis_japonica.C_japonica-7.0.1.22.dna_sm.stable.fa',
    G_gallus = '/cbio/grlab/share/databases/genomes/G_gallus/ensembl_release-69/Gallus_gallus.WASHUC2.69_stable.fa.bz2',
    M_musculus = '/cbio/grlab/share/databases/genomes/M_musculus/ensembl_release-69/Mus_musculus.GRCm38.69_stable.fa.bz2',  
    V_vinifera = '/cbio/grlab/share/databases/genomes/V_vinifera/STARgenome/Vvinifera_145.fa',
    A_mellifera = '/cbio/grlab/share/databases/genomes/A_mellifera/ensembl_release-20/apiMel3_ucsc_dnaMaskchrom.fasta',
    B_taurus = '/cbio/grlab/share/databases/genomes/B_taurus/STARgenome/ensembl_release-69.fa',
    C_rubella = '/cbio/grlab/share/databases/genomes/C_rubella/phytozome_v9.0/Crubella_183.fa.gz',
    D_rerio = '/cbio/grlab/share/databases/genomes/D_rerio/ensembl_release-69/Danio_rerio.Zv9.69.dna.toplevel.fa.bz2',
    G_max = '/cbio/grlab/share/databases/genomes/G_max/phytozome_v9.0/Gmax_189_filter.fa',
    M_truncatula = '/cbio/grlab/share/databases/genomes/M_truncatula/STARgenome/Mtruncatula_198.fa',
    P_pacificus = '/cbio/grlab/share/databases/genomes/P_pacificus/STARgenome/Pristionchus_pacificus.P_pacificus-5.0.22.dna_sm.stable.fa',
    S_scrofa = '/cbio/grlab/share/databases/genomes/S_scrofa/ensembl_release-69/Sus_scrofa.Sscrofa10.2.69.dna_filter.fasta',
    X_tropicalis = '/cbio/grlab/share/databases/genomes/X_tropicalis/ensembl_release-69/Xenopus_tropicalis.JGI_4.2.69.dna_stable.fasta.bz2',
    C_sativus = '/cbio/grlab/share/databases/genomes/C_sativus/phytozome_v9.0/Csativus_122_filtered.fa',
    D_simulans = '/cbio/grlab/share/databases/genomes/D_simulans/ensembl_release-22/Drosophila_simulans.WUGSC1.22.dna_sm.dna.fa',
    H_sapiens = '/cbio/grlab/share/databases/genomes/H_sapiens/hg19_bowtie/hg19.fa',
    O_anatinus = '/cbio/grlab/share/databases/genomes/O_anatinus/ensembl_release-69/Ornithorhynchus_anatinus.OANA5.69-filtered_dna.fa',
    P_troglodytes = '/cbio/grlab/share/databases/genomes/P_troglodytes/STARgenome/Pan_troglodytes.CHIMP2.1.4.69_stable.fa',
    S_tuberosum = '/cbio/grlab/share/databases/genomes/S_tuberosum/phytozome_v9.0/Stuberosum_206.fa',
    Z_mays = '/cbio/grlab/share/databases/genomes/Z_mays/phytozome_v9.0/Zmays_181.fa',
    A_thaliana = '/cbio/grlab/share/databases/genomes/A_thaliana/arabidopsis_tair10/sequences/TAIR9_chr_all.fas',
    C_elegans = '/cbio/grlab/share/databases/genomes/C_elegans/ensembl_release-69/Caenorhabditis_elegans.WBcel215.69.dna.toplevel.fa',
    D_discoideum = '/cbio/grlab/share/databases/genomes/D_discoideum/STARgenome/Dictyostelium_discoideum.dictybase.01.21.dna.toplevel.fa',
    D_yakuba = '/cbio/grlab/share/databases/genomes/D_yakuba/STARgenome/Drosophila_yakuba.dyak_r1.3_FB2008_07.22.dna_sm_stable.fasta',
    O_latipes = '/cbio/grlab/share/databases/genomes/O_latipes/STARgenome/Oryzias_latipes.MEDAKA1.74.dna_rm.toplevel.fa',
    R_norvegicus = '/cbio/grlab/share/databases/genomes/R_norvegicus/ensembl_release-69/Rattus_norvegicus.RGSC3.4.69.dna.toplevel.fa.bz2',
    C_briggsae = '/cbio/grlab/share/databases/genomes/C_briggsae/STARgenome/Caenorhabditis_briggsae.CB4.22.dna_sm_stable.fasta',
    C_brenneri = '/cbio/grlab/share/databases/genomes/C_brenneri/STARgenome/Caenorhabditis_brenneri.C_brenneri-6.0.1b.22.dna_sm_stable.fa',
    C_remanei = '/cbio/grlab/share/databases/genomes/C_remanei/STARgenome/Caenorhabditis_remanei.C_remanei-15.0.1.22.dna_sm_stable.fa',
    D_pseudoobscura = '/cbio/grlab/share/databases/genomes/D_pseudoobscura/STARgenome/Drosophila_pseudoobscura.HGSC2.22.dna_sm_stable.fasta',
    T_nigroviridis = '/cbio/grlab/share/databases/genomes/T_nigroviridis/ensembl_release-69/Tetraodon_nigroviridis.TETRAODON8.69.dna.toplevel.fa'
    )

base_path_gff = '/cbio/grlab/nobackup/SignalPrediction/SRA-rnaseq/%s/trans_pred/filtered_genes.gff' % org_name

generate_genome_seq_labels.main(org_details[org_name], base_path_gff, 'tss', 6300, 5000, 15000, 1200)

#generate_genome_seq_labels.main(org_details[org_name], base_path_gff, 'cleave', 1200)
#generate_genome_seq_labels.main(org_details[org_name], base_path_gff, 'cdsstop', 1600)
#generate_genome_seq_labels.main(org_details[org_name], base_path_gff, 'tss', 1200)
#generate_genome_seq_labels.main(org_details[org_name], base_path_gff, 'tis', 2000)
