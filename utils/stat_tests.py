#!/usr/bin/env python 
"""
statistical test for determining the performance of learning algorithms
"""

import numpy
import pandas 
from scipy import stats 
from collections import defaultdict

from utils import compressed_pickle


def best_test_perf_with_eval(xv_result_file, methods=None, org_names=None):
    """
    better test performance value vased on the model validation performance for each XV.
    """
    
    data = compressed_pickle.load(xv_result_file) 

    if methods is None:
        methods = data.keys()

    if org_names is None:
        org_names = data[methods[0]][0].keys()
    
    best_test_perf = defaultdict()

    for m_idx, m in enumerate(methods):
        for n_idx, n in enumerate(org_names):

            assert data[m][0][n].shape == data[m][1][n].shape

            num_splits, num_params = data[m][0][n].shape
            best_eval_param_idx_for_each_split = numpy.argmax(data[m][0][n], axis=1) ##top val for each XV
            
            assert len(best_eval_param_idx_for_each_split) == num_splits

            tmp_test = []
            for s_idx, p_idx in enumerate(best_eval_param_idx_for_each_split):
                tmp_test.append(data[m][1][n][s_idx][p_idx])

            assert len(tmp_test) == num_splits

            if n in best_test_perf:
                best_test_perf[n][m] = numpy.array(tmp_test)
            else:
                best_test_perf[n] = {m:numpy.array(tmp_test)}

    return best_test_perf


def compute_friedman_test(best_test_score):
    """
    """
   
    ## test pair combinations 
    pairs_test = [('individual', 'union', 'mtl'), ('individual', 'union', 'mtmkl'), \
        ('individual', 'mtl', 'mtmkl'), ('union', 'mtl', 'mtmkl')]
    org_names = best_test_score.keys()

    ##ttest_p_val = numpy.zeros((len(org_names), len(pairs_test)))
    ttest_p_val = numpy.zeros((len(pairs_test), len(org_names)))

    for org_idx, org_code in enumerate(org_names):
        meth_perf = best_test_score[org_code]

        for pair_idx, rel_pair in enumerate(pairs_test):
            t_stats, p_val = stats.friedmanchisquare(meth_perf[rel_pair[0]], meth_perf[rel_pair[1]], meth_perf[rel_pair[2]])
            
            ##ttest_p_val[org_idx, pair_idx] = p_val
            ttest_p_val[pair_idx, org_idx] = p_val
        
    
    ##df_pval = pandas.DataFrame(ttest_p_val, columns=pairs_test, index=org_names)
    df_pval = pandas.DataFrame(ttest_p_val, columns=org_names, index=pairs_test)
    
    return df_pval 


def compute_wilcoxon_test(best_test_score):
    """
    """
    
    pairs_test = [('individual', 'union'), ('individual', 'mtl'), \
        ('individual', 'mtmkl'), ('union', 'mtl'), ('union', 'mtmkl'), \
        ('mtl', 'mtmkl')]
    org_names = best_test_score.keys()

    ##ttest_p_val = numpy.zeros((len(org_names), len(pairs_test)))
    ttest_p_val = numpy.zeros((len(pairs_test), len(org_names)))

    for org_idx, org_code in enumerate(org_names):
        meth_perf = best_test_score[org_code]

        for pair_idx, rel_pair in enumerate(pairs_test):
            t_stats, p_val = stats.wilcoxon(meth_perf[rel_pair[0]], meth_perf[rel_pair[1]])
            
            ##ttest_p_val[org_idx, pair_idx] = p_val
            ttest_p_val[pair_idx, org_idx] = p_val
        
    
    ##df_pval = pandas.DataFrame(ttest_p_val, columns=pairs_test, index=org_names)
    df_pval = pandas.DataFrame(ttest_p_val, columns=org_names, index=pairs_test)
    
    return df_pval 


def compute_paired_ttest(best_test_score):
    """
    """
    
    pairs_test = [('individual', 'union'), ('individual', 'mtl'), \
        ('individual', 'mtmkl'), ('union', 'mtl'), ('union', 'mtmkl'), \
        ('mtl', 'mtmkl')]
    org_names = best_test_score.keys()

    ttest_p_val = numpy.zeros((len(org_names), len(pairs_test)))
    #ttest_p_val = numpy.zeros((len(pairs_test), len(org_names)))

    for org_idx, org_code in enumerate(org_names):
        meth_perf = best_test_score[org_code]

        for pair_idx, rel_pair in enumerate(pairs_test):
            t_stats, p_val = stats.ttest_rel(meth_perf[rel_pair[0]], meth_perf[rel_pair[1]])
            
            ttest_p_val[org_idx, pair_idx] = p_val
            #ttest_p_val[pair_idx, org_idx] = p_val
        
    
    df_pval = pandas.DataFrame(ttest_p_val, columns=pairs_test, index=org_names)
    #df_pval = pandas.DataFrame(ttest_p_val, columns=org_names, index=pairs_test)
    
    return df_pval 
