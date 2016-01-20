#!/usr/bin/env python 
"""
training a model with single or multiple organisms 
"""

## standard modules 
import time 
import uuid
import numpy
from Bio import SeqIO

## custom modules
import helper 
from utils import compressed_pickle
from promoter_kernel import ShogunPredictor


def load_examples_from_fasta(signal, org, data_path):
    """
    load examples from fasta file
    
    @args signal: genomic signal type (default: tss) 
    @type signal: str 
    @args org: organism name (ex: A_thaliana) 
    @type org: str 
    @args data_path: file path for training data points 
    @type data_path: str 
    """
    
    ## defining the data point complete path   
    fn_pos = "%s/%s_sig_%s_example.fa" % (data_path, signal, "pos")
    fn_neg = "%s/%s_sig_%s_example.fa" % (data_path, signal, "neg")
    print "loading: \n %s \n %s" % (fn_pos, fn_neg) 

    ## parse examples from fasta file
    xt_pos = [str(rec.seq) for rec in SeqIO.parse(fn_pos, "fasta")]
    xt_neg = [str(rec.seq) for rec in SeqIO.parse(fn_neg, "fasta")]

    labels = [+1] * len(xt_pos) + [-1] * len(xt_neg)
    examples = xt_pos + xt_neg

    ## some log information 
    print("organism: %s, signal %s,\t num_labels: %i,\t num_examples %i,\t num_positives: %i,\t num_negatives: %i" %  (org, signal, len(labels), len(examples), len(xt_pos), len(xt_neg)))

    examples_shuffled, labels_shuffled = helper.coshuffle(examples, labels)
    ret = {"examples": numpy.array(examples_shuffled), "labels": numpy.array(labels_shuffled)}

    return ret


def train_wdspeck_svm(org_code, signal="tss", data_path="SRA-rnaseq"):
    """
    train SVM based on the examples from different sources 

    @args org_code: organism name (ex: A_thaliana)
    @type org_code: str 
    @args signal: genomic signal type (default: tss) 
    @type signal: str 
    @args data_path: file path for training data points 
    @type data_path: str 
    """
    t0 = time.time()

    ## loading data
    data = load_examples_from_fasta(signal, org_code, data_path)
    assert(len(data["examples"]) == len(data["labels"]))

    ## split the data 
    train_examples = data["examples"]
    train_labels = data["labels"]

    ## set parameters TODO 
    param = {}
    param["cost"] = 1.0
    param["degree"] = 4 
    param["degree_spectrum"] = 4
    param["center_pos"] = 1200
    param["center_offset"] = 50 
    param["shifts"] = 32
    param["kernel_cache"] = 10000

    ## invoke training
    svm = ShogunPredictor(param)
    svm.train(train_examples, train_labels)

    ## save the model 
    fname = "%s_%s_model_%s" % (org_code, signal, uuid.uuid1()) 
    compressed_pickle.save(fname, svm) 
    print("saving the model in file %s" % fname)

    time_taken = time.time() - t0
    print("time taken for the experiment: ", time_taken)

    return fname 


def get_orgdb(fn):
    """
    retrive organism name from txt file 
    """

    org_names = []
    for line in open(fn):
        if line.strip() != "" and not line.startswith("#"):
            org_names.append(line.strip())

    org_names.sort()
    print org_names

    return org_names


def train_combined_wdspeck_svm(org_list_file, signal="tss", data_path="SRA-seq"):
    """
    train a global classifier for multiple organisms

    @args org_list_file: organism name in a text file 
    @type org_list_file: str 
    @args signal: genomic signal type (default: tss) 
    @type signal: str 
    @args data_path: file path for training data points 
    @type data_path: str 
    """

    t0 = time.time()
    ## getting multiple oragnisms 
    ORG_LIST = get_orgdb(org_list_file)

    ## loading data
    train_examples = [] 
    train_labels = []
    for ORG_NAME in ORG_LIST: 
        local_data_path = "%s/%s/set_1" % (data_path, ORG_NAME) ## FIXME common data path
        data = load_examples_from_fasta(signal, ORG_NAME, local_data_path)
        assert(len(data["examples"]) == len(data["labels"]))

        ## split the data 
        train_examples.extend(data["examples"])
        train_labels.extend(data["labels"])

    ## set parameters TODO
    param = {}
    param["cost"] = 1.0
    param["degree"] = 4 
    param["degree_spectrum"] = 4
    param["center_pos"] = 1200
    param["center_offset"] = 50 
    param["shifts"] = 32
    param["kernel_cache"] = 10000

    ## invoke training
    svm = ShogunPredictor(param)
    svm.train(train_examples, train_labels)

    ## save the model 
    fname = "%s_model_%s" % (signal, uuid.uuid1()) 
    compressed_pickle.save(fname, svm) 
    print("saving the model in file %s" % fname)

    time_taken = time.time() - t0
    print("time taken for the experiment: ", time_taken)

    return fname 
