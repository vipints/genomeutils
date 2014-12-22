#!/usr/bin/env python 
"""
Understand the genome quality like contig numbers and non nucleotide
characters. 

Usage: 
Requirement: 
"""

import sys, os, re 
import gzip 
import bz2
from Bio import SeqIO

def _open_file(fname):
    """returns the file handler
    """
    if os.path.exists(fname): 
        if os.path.splitext(fname)[1] == ".gz":
            FH = gzip.open(fname, 'rb')
        elif os.path.splitext(fname)[1] == ".bz2":
            FH = bz2.BZ2File(fname, 'rb')
        else:
            FH = open(fname, 'rU')
    else:
        sys.exit('File not found ' + fname) 
    return FH 

def seq_analysis(fhd):
    """read the genome sequence 
    """
    seq_info = dict()
    for rec in SeqIO.parse(fhd, "fasta"):
        seq_info[rec.id] = len(rec.seq)
    fhd.close()

def find_org_files(base_path, genus):
    """
    get the genome sequence file and get the details 
    """
    base_path = os.path.realpath(base_path)
    spec_file_path = base_path + '/' + genus +'/ensembl_release-69/'
    dir_contents = os.listdir(spec_file_path)

    #ensembl dir contains fasta dna.toplevel.fa and gtf 
    fastafn, gtfname = None, None 
    for filename in dir_contents:
        if re.search(r'dna.toplevel', filename):
            fastafn = filename 
        if re.search(r'.gtf', filename):
            gtfname = filename 
    
    #check the genome sequnce to get the quality of sequence
    fh = _open_file(spec_file_path + fastafn)
    seq_analysis(fh)

def init_org():
    org = dict(
            A_melanoleuca = 1, 
            A_carolinensis  = 1,
            B_taurus = 1,
            C_elegans = 1, 
            C_jacchus = 1,
            C_familiaris = 1,
            C_porcellus = 1,
            C_hoffmanni = 1,
            C_intestinalis = 1,
            C_savignyi = 1,
            D_rerio = 1,
            D_novemcinctus = 1,
            D_ordii = 1,
            D_melanogaster = 1,
            E_telfairi = 1,
            E_caballus = 1,
            E_europaeus = 1,
            F_catus = 0,
            G_morhua = 1,
            G_gallus = 0,
            G_aculeatus = 0,
            G_gorilla = 0,
            H_sapiens = 1,
            I_tridecemlineatus = 1,
            L_chalumnae = 0,
            L_africana = 0,
            M_mulatta = 0,
            M_eugenii = 1,
            M_gallopavo = 0,
            M_murinus = 0,
            M_domestica = 0,
            M_musculus = 1,
            M_putorius_furo = 1,
            M_lucifugus = 1,
            N_leucogenys = 1,
            O_princeps = 1,
            O_niloticus = 1,
            O_anatinus = 1,
            O_cuniculus = 1,
            O_latipes = 1,
            O_garnettii = 1,
            P_troglodytes = 1,
            P_sinensis = 1,
            P_marinus = 1,
            P_abelii = 1,
            P_capensis = 1,
            P_vampyrus = 1,
            R_norvegicus = 1,
            S_cerevisiae = 1,
            S_harrisii = 1,
            S_araneus = 1,
            S_scrofa = 1,
            T_guttata = 1,
            T_rubripes = 1,
            T_syrichta = 1,
            T_nigroviridis = 1,
            T_belangeri = 1,
            T_truncatus = 1,
            V_pacos = 1,
            X_tropicalis = 1,
            X_maculatus = 1
            )
    return org 

if __name__=="__main__":
    try:
        base_dir = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    
    species = init_org()
    for spname in sorted(species):
        print spname 
        find_org_files(base_dir, spname)
        break
