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


def get_trsk_genes(gtf_file):
    """
    parse the trsk predictions 
    create the df with length 
    only select the genes with CDS region 
    """
    
    anno_db = GFFParser.Parse(gtf_file) 

    #cnt = 0 
    #for feat in anno_db:
    #    if not feat['cds_exons'][0].any(): ## trsk produce single transcript 
    #        cnt += 1 
    #genes_w_cds = len(anno_db) - cnt 

    trans_len = np.zeros(len(anno_db, 2))
    genes = [] 

    for idx, feat in enumerate(anno_db):
        cds_len = 0 
        for exc in feat['cds_exons'][0]:
            cds_len += exc[1]-exc[0]

        trans_len[idx, 0] = feat['stop']-feat['start']
        trans_len[idx, 1] = cds_len
        genes.append(feat['name'])
        
    df_len_dis_genes = pd.DataFrame(trans_len, columns=['gene_len', 'cds_len'], index=genes)


if __name__=="__main__":
    fname = "../../SRA-rnaseq/H_sapiens/trans_pred/H_sapiens_trsk_genes.gff"
    get_trsk_genes(fname)
