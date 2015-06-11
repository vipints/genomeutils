#!/usr/bin/env python 
"""
translate the CDS region predicted by trsk program 
"""

import os 
import sys 
import numpy as np 
import pandas as pd 
from Bio import SeqIO 
from gfftools import GFFParser, helper

def translate_trsk_genes(gtf_file, fas_file):
    """
    translate the trsk genes to protein sequence 
    """
    #from Bio.Seq import Seq
    #from Bio.Alphabet import generic_rna, generic_dna

    anno_db = GFFParser.Parse(gtf_file) 
    total_genes = len(anno_db) 

    cds_idx = [] # deleting the empty cds lines  
    for idp, feat in enumerate(anno_db):
        if not feat['cds_exons'][0].any():
            cds_idx.append(idp) 
    anno_db = np.delete(anno_db, cds_idx) 
    genes_with_cds = len(anno_db) 

    ## genome sequence file reading 
    fasFH = helper.open_file(fas_file) 

    for rec in SeqIO.parse(fasFH, "fasta"):
        for idx, feature in enumerate(anno_db):
            if rec.id == feat['chr']:
                ## iterate over cds_exons
                cds_seq = ''
                for ex in feature['cds_exons'][0]:## single transcript by TSkim 
                    cds_seq += rec.seq[ex[0]-1:ex[1]]
                
                if feature['strand'] == '-':
                    cds_seq = cds_seq.reverse_complement()

                sys.stdout.write(str(cds_seq.translate()) + "\n")

        # FIXME need an efficient way to translate multiple gene 
        # iterate over chromosome
    fasFH.close()


def trsk_gene_len_dist(gtf_file):
    """
    plotting the histograms bases on the genes and CDS length
    """
    import matplotlib.pyplot as plt 

    anno_db = GFFParser.Parse(gtf_file) 

    cds_idx = [] # deleting the empty cds lines  
    for idp, feat in enumerate(anno_db):
        if not feat['cds_exons'][0].any():
            cds_idx.append(idp) 

    anno_db = np.delete(anno_db, cds_idx) 
    
    trans_len = np.zeros((len(anno_db), 2))
    genes = [] 

    for idx, feat in enumerate(anno_db):
        cds_len = 0 
        for exc in feat['cds_exons'][0]:
            cds_len += exc[1]-exc[0]

        trans_len[idx, 0] = feat['stop']-feat['start']
        trans_len[idx, 1] = cds_len
        genes.append(feat['name'])
    
    ## gene, cds length information 
    df_len_dis_genes = pd.DataFrame(trans_len, columns=['gene_len', 'cds_len'], index=genes)

    ## plotting the gene length based on the bins of gene length  
    gene_length = trans_len[:,0] ## gene length from the matrix 

    freq, bins = np.histogram(gene_length, bins=10, range=None, normed=False, weights=None)
    bins = np.delete(bins, 10) 

    df_gene_len_bin = pd.DataFrame(freq, columns=['gene_frequency'], index=bins) 
    plt.figure() 
    df_gene_len_bin.plot(kind="bar")
    #plt.savefig()

    ## plotting the cds length distribution
    cds_length = trans_len[:,1] ## cds length distribution 
    freq, bins = np.histogram(cds_length, bins=10, range=None, normed=False, weights=None)
    bins = np.delete(bins, 10) 
    df_cds_len_bin = pd.DataFrame(freq, columns=['cds_frequency'], index=bins) 
    plt.figure() 
    df_cds_len_bin.plot(kind="bar")
    #plt.savefig("hist_cds_len.pdf") 


if __name__=="__main__":
    gname = "hs_chr22.gff"
    fas = "hg19.fa"
    translate_trsk_genes(gname, fas) 
