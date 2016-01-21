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




def predict_around_region(args_list):
    """
    predict around the center offset region of the signal
    """

    start_scan, stop_scan, model, data_set = args_list 
    trimed_seq_list = [] 

    for datum in data_set:
        tss_score = numpy.zeros(100) ## predefined flanking region of 100 nts 

        ## predicting the near regions of the signal, here it is -/+ 50 nucleotides
        for idy, shifted_cen in enumerate(xrange(start_scan, stop_scan)):
            model.param["center_pos"] = shifted_cen
            out = model.predict(datum)
            tss_score[idy] = out[0] 

        trimed_seq_list.append(tss_score)
    
    return trimed_seq_list


def reduce_pred_score(flat_result_list):
    """
    collect the results from different worker
    """
    
    idx = 0 
    merged_arr = numpy.zeros((100, 100))

    for flat_result in flat_result_list:
        for element in flat_result:
            merged_arr[idx] = element
            idx += 1 

    return merged_arr


def load_model(svm_file):
    """
    load the model from a pickled file
    """

    import bz2 
    import cPickle 

    fh = bz2.BZ2File(svm_file, "rb")
    svm = cPickle.load(fh) 
    fh.close() 

    return svm


def load_examples_from_fasta(signal, ex_type, org, data_path):
    """
    load examples from fasta file
    
    @args signal: genomic signal type (default: tss) 
    @type signal: str 
    @args ex_type: example type (default: pos) 
    @type ex_type: str 
    @args org: organism name (ex: A_thaliana) 
    @type org: str 
    @args data_path: file path for training data points 
    @type data_path: str 
    """
    
    fn_pos = "%s/%s_sig_%s_example.fa" % (data_path, signal, ex_type)
    print("loading: \n %s" % (fn_pos))

    # parse file
    xt_pos = [str(rec.seq) for rec in SeqIO.parse(fn_pos, "fasta")]
    labels =  [+1] * len(xt_pos)
    examples =  xt_pos

    print("organism: %s, signal %s,\t num_labels: %i,\t num_examples %i,\t num_positives: %i" %  (org, signal, len(labels), len(examples), len(xt_pos)))

    examples_shuffled, labels_shuffled = helper.coshuffle(examples, labels)
    ret = {"examples": numpy.array(examples_shuffled), "labels": numpy.array(labels_shuffled)}

    return ret


def recenter_examples(args_list):
    """
    According to the local maximum prediction output score, 
    manually fix each example sequence signal site position.
    """

    start_scan, stop_scan, model, data_set = args_list 

    trimed_seq_list = [] 
    cen_flank_region = stop_scan - start_scan
    site_position = model.param["center_pos"]

    for datum in data_set:## wrapping multiple sequence 
        tss_score = numpy.zeros(cen_flank_region) 

        ## predicting the near regions of the signal - center site flanking region 
        for idy, shifted_cen in enumerate(xrange(start_scan, stop_scan)):
            model.param["center_pos"] = shifted_cen
            out = model.predict(datum)
            tss_score[idy] = out[0] 

        max_pred_out = numpy.argmax(tss_score)
        left_offset = site_position - start_scan
        
        if max_pred_out < left_offset: ## center_offset left 
            rel_cen_pos = (site_position - left_offset) + max_pred_out ## center_pos 
        else:
            rel_cen_pos = site_position + (max_pred_out-left_offset)
        
        trim_reg_site = site_position - cen_flank_region ## 1200 - 100 = 1100
        trimed_seq = datum[0][rel_cen_pos-trim_reg_site:rel_cen_pos+trim_reg_site] 
        trimed_seq_list.append(trimed_seq)

    return trimed_seq_list



def shift_signal_position(svm_file, org, example_type="pos", signal="tss", data_path="SRA-rnaseq"):
    """
    manually look at the position around the original position 
    """

    ## loading data
    data = load_examples_from_fasta(signal, example_type, org, data_path)
    assert(len(data["examples"]) == len(data["labels"]))

    ## unpack the model
    model = load_model(svm_file)
   
    ## getting the model information 
    center_pos = model.param["center_pos"]
    center_offset = model.param["center_offset"]
    
    print("model - center pos: %i, center reg: %i" % (center_pos, center_offset))
    
    ## the recentering the center regions 
    start_scan = center_pos-center_offset
    stop_scan = center_pos+center_offset

    cnt = 0
    data_set = [] 
    argument_list = [] 

    ## get the individual examples to recenter the signal position manually
    for idx, single_example in enumerate(data["examples"]):  
        datum = [single_example]
        cnt += 1

        if cnt % 10 == 0: ## packing 10 seq to one job 
            data_set.append(datum)

            arg = [start_scan, stop_scan, model, data_set]
            argument_list.append(arg)

            data_set = [] 
        else:
            data_set.append(datum)
    
    local = False ## switch between local and compute cluster 
    ## cluster compute options   
    cluster_resource = {'pvmem':'8gb', 'pmem':'8gb', 'mem':'8gb', 'vmem':'8gb','ppn':'1', 'nodes':'1', 'walltime':'24:00:00'}

    ## job dispatching 
    intm_ret = pg.pg_map(recenter_examples, argument_list, param=cluster_resource, local=local, maxNumThreads=1, mem="8gb") 

    if task_type:
        print "Done with computation"

        fixed_example_seq = reduce_modified_seq(intm_ret) 
        print "Done reducing the results"

        write_fasta_rec(fixed_example_seq, signal) 

    else:
        intm_ret = pg.pg_map(predict_around_region, argument_list, param=cluster_resource, local=local, maxNumThreads=2, mem="4gb") 
        print "Done with computation"

        pred_out_val = reduce_pred_score(intm_ret) 
        print "Done reducing the results"
    
        ## save the scores 
        fname = "%s_pred_score_%s" % (signal, uuid.uuid1()) 
        compressed_pickle.save(fname, pred_out_val) 

        print( "saving the score in file %s" % fname )


def main():

    org_code = "H_sapiens"


if __name__ == "__main__":
    main()

