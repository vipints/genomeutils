#!/usr/bin/env python
"""
module to recenter the labels based on the manual shifting of tss signal region.
"""

import uuid
import numpy
from Bio import SeqIO

from data_processing import helper
from utils import compressed_pickle

from promoter_kernel import ShogunPredictor

import libpyjobrunner as pg


def load_data_from_fasta(signal, org, data_path):
    """
    load examples and labels from fasta file
    """
    
    #fn_pos = "%s/%s_sig_%s_example_900.fa" % (data_path, signal, "pos")
    fn_neg = "%s/%s_sig_%s_example_1800.fa" % (data_path, signal, "neg")

    #print "loading: \n %s \n %s" % (fn_pos, fn_neg) 
    #print "loading: \n %s" % (fn_pos) 
    print "loading: \n %s" % (fn_neg) 

    # parse file
    #xt_pos = [str(rec.seq) for rec in SeqIO.parse(fn_pos, "fasta")]
    xt_neg = [str(rec.seq) for rec in SeqIO.parse(fn_neg, "fasta")]

    #labels = [+1] * len(xt_pos) + [-1] * len(xt_neg)
    #labels =  [+1] * len(xt_pos)
    labels =  [-1] * len(xt_neg)

    #examples = xt_pos + xt_neg
    #examples =  xt_pos
    examples =  xt_neg

    #print("organism: %s, signal %s,\t num_labels: %i,\t num_examples %i,\t num_positives: %i,\t num_negatives: %i" %  (org, signal, len(labels), len(examples), len(xt_pos), len(xt_neg)))
    #print("organism: %s, signal %s,\t num_labels: %i,\t num_examples %i,\t num_positives: %i" %  (org, signal, len(labels), len(examples), len(xt_pos)))
    print("organism: %s, signal %s,\t num_labels: %i,\t num_examples %i,\t num_negatives: %i" %  (org, signal, len(labels), len(examples), len(xt_neg)))

    examples_shuffled, labels_shuffled = helper.coshuffle(examples, labels)

    ret = {"examples": numpy.array(examples_shuffled), "labels": numpy.array(labels_shuffled)}

    return ret


def reduce_modified_seq(flat_result):
    """ return the result 
    """
    
    return flat_result


def write_fasta_rec(seq_list_total, signal):
    """ write fasta file from the seq records 
    """
    seq_type = "+1" 
    fname = "trimmed_%s_pos_examples.fa" % signal

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    out_fh = open(fname, "w")
    for seq_list in seq_list_total: 
        for seq_rec in seq_list:
            faseq_out = SeqRecord(Seq(seq_rec), id="example_%s" % uuid.uuid1(), description=seq_type)
            out_fh.write(faseq_out.format('fasta'))
    out_fh.close() 


def load_examples_from_fasta(signal, org, data_path):
    """
    load examples from fasta file
    """
    
    fn_pos = "%s/%s_sig_%s_example.fa" % (data_path, signal, "pos")
    fn_neg = "%s/%s_sig_%s_example.fa" % (data_path, signal, "neg")
    print "loading: \n %s \n %s" % (fn_pos, fn_neg) 

    # parse file
    xt_pos = [str(rec.seq) for rec in SeqIO.parse(fn_pos, "fasta")]
    xt_neg = [str(rec.seq) for rec in SeqIO.parse(fn_neg, "fasta")]

    labels = [+1] * len(xt_pos) + [-1] * len(xt_neg)
    examples = xt_pos + xt_neg

    print("organism: %s, signal %s,\t num_labels: %i,\t num_examples %i,\t num_positives: %i,\t num_negatives: %i" %  (org, signal, len(labels), len(examples), len(xt_pos), len(xt_neg)))

    examples_shuffled, labels_shuffled = helper.coshuffle(examples, labels)
    ret = {"examples": numpy.array(examples_shuffled), "labels": numpy.array(labels_shuffled)}

    return ret


def train_shifted_wdk_svm(org_code, signal="tss", data_path="SRA-rnaseq"):
    """
    train SVM based on the examples from different sources 
    """

    import time 
    t0 = time.time()

    ## loading data
    data = load_examples_from_fasta(signal, org_code, data_path)
    assert(len(data["examples"]) == len(data["labels"]))

    ## split the data 
    train_examples = data["examples"]
    train_labels = data["labels"]

    ## set parameters
    param = {}
    param["cost"] = 1.0
    param["degree"] = 8
    param["degree_spectrum"] = 4
    param["center_pos"] = 1200
    param["center_offset"] = 50 
    param["shifts"] = 32
    param["kernel_cache"] = 8000

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


def plot_tss_score(tss_score_file):
    """
    """

    import bz2 
    import cPickle 
    import matplotlib.pyplot as plt

    ## local maxima and local minima 
    from scipy.signal import argrelextrema

    fh = bz2.BZ2File(tss_score_file, "rb")

    tss_score = cPickle.load(fh)
    fh.close()

    #argrelextrema(tss_score[0], numpy.greater)
    #plt.annotate('local max', xy=(2, 1), xytext=(3, 1.5),
    #            arrowprops=dict(facecolor='black', shrink=0.05),
    #                        )

    dt = 1.0
    t = numpy.arange(0, len(tss_score[0]), dt)

    ## plotting two lines
    plt.plot(t, tss_score[0], 'b-', t, tss_score[1], 'g-')
    plt.xlim(0.0,100.0)

    plt.xlabel('-/+ 50 nucleotides flanking region of tss signal')
    plt.ylabel('model prediction score')
    plt.grid(True)

    #plt.show()
    out_file = "arts.pdf"
    plt.savefig(out_file)


def predict_and_recenter(args_list):
    """
    fix the position
    recenter the signal position according to the local max pred out 
    """

    start_scan, stop_scan, model, data_set = args_list 
    trimed_seq_list = [] 
    
    for datum in data_set:## wrapping multiple sequence 

        tss_score = numpy.zeros(100) ## predefined flanking region of 100 nts 

        ## predicting the near regions of the signal, here it is -/+ 50 nucleotides
        for idy, shifted_cen in enumerate(xrange(start_scan, stop_scan)):
            model.param["center_pos"] = shifted_cen
            out = model.predict(datum)
            tss_score[idy] = out[0] 

        max_pred_out = numpy.argmax(tss_score)
        center_pos = 1200
        left_offset = center_pos - start_scan
        
        if max_pred_out < left_offset: ## center_offset left 
            rel_cen_pos = (center_pos - left_offset) + max_pred_out ## center_pos 
        else:
            rel_cen_pos = center_pos + (max_pred_out-left_offset)

        trimed_seq = datum[0][rel_cen_pos-1100:rel_cen_pos] + datum[0][rel_cen_pos:(rel_cen_pos+1100)] 
        trimed_seq_list.append(trimed_seq)

    return trimed_seq_list


def manual_pos_shift(svm_file, org, signal="tss", data_path="SRA-rnaseq"):
    """
    manually look at the position around the original position 
    """

    ## loading data
    data = load_examples_from_fasta(signal, org, data_path)
    assert(len(data["examples"]) == len(data["labels"]))

    ## unpack the model
    import bz2 
    import cPickle 

    fh = bz2.BZ2File(svm_file, "rb")
    model = cPickle.load(fh) 
    fh.close() 
    
    ## getting the model information 
    center_pos = model.param["center_pos"]
    center_offset = model.param["center_offset"]
    
    print("model - center pos: %i, center reg: %i" % (center_pos, center_offset))

    start_scan = center_pos-center_offset
    stop_scan = center_pos+center_offset

    cnt = 0
    data_set = [] 
    argument_list = [] 

    label_type = -1 ## label_type will be +1/-1

    ## get the individual examples to recenter the signal position manually
    for idx, single_example in enumerate(data["examples"]):  

        datum = [single_example]
        label_info = data['labels'][idx]
        
        if label_info != label_type: 
            cnt += 1

            if cnt % 10 == 0: ## packing 10 seq to one job 
                data_set.append(datum)

                arg = [start_scan, stop_scan, model, data_set]
                argument_list.append(arg)

                data_set = [] 
            else:
                data_set.append(datum)
                
    
    local = False 
    cluster_resource = {'pvmem':'4gb', 'pmem':'4gb', 'mem':'4gb', 'vmem':'4gb','ppn':'1', 'nodes':'1', 'walltime':'4:00:00'}
    task_type = 0 # 1 recenter seq, 0 predict score  

    if task_type:

        intm_ret = pg.pg_map(predict_and_recenter, argument_list, param=cluster_resource, local=local, maxNumThreads=2, mem="4gb") 
        print "Done with computation"

        fixed_example_seq = reduce_modified_seq(intm_ret) 
        print "Done reducing the results"

        write_fasta_rec(fixed_example_seq, signal) 

    else:
        
        intm_ret = pg.pg_map(predict_around_region, argument_list, param=cluster_resource, local=local, maxNumThreads=2, mem="4gb") 
        print "Done with computation"

        #import ipdb 
        #ipdb.set_trace() 

        pred_out_val = reduce_result(intm_ret) 
        print "Done reducing the results"
    
        ## save the scores 
        fname = "%s_pos_centeroff2K_score_%s" % (signal, uuid.uuid1()) 
        compressed_pickle.save(fname, pred_out_val) 
        print( "saving the score in file %s" % fname )



def reduce_result(flat_result):
    """
    collect the results from different worker
    """

    merged_arr = numpy.zeros((len(flat_result), len(flat_result[0])))
    
    for idx, element in enumerate(flat_result):
        merged_arr[idx] = element

    return merged_arr


def main():

    org_code = "H_sapiens"
    model_file_name = train_shifted_wdk_svm(org_code)

    model_file_name =  "tss_28a9ce8c-528f-11e5-b86a-90e2ba3a745c"
    manual_pos_shift(model_file_name, org_code)
 

if __name__ == "__main__":
    main()

