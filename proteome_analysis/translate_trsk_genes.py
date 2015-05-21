#!/usr/bin/env python 
"""
translate the CDS region predicted by trsk program 
"""

import os 
import sys 
import numpy as np 
import pandas as pd 
from Bio import SeqIO 
from gfftools import GFFParser
import matplotlib.pyplot as plt 


def trsk_gene_len_dist(gtf_file):
    """

    """
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

    df_gene_len_bin = pd.DataFrame(freq, columns=['frequency'], index=bins) 
    plt.figure() 

    df_gene_len_bin.plot(kind="bar")
    #plt.savefig()


if __name__=="__main__":
    fname = "../../SRA-rnaseq/H_sapiens/trans_pred/H_sapiens_trsk_genes.gff"
    trsk_gene_len_dist(fname)
