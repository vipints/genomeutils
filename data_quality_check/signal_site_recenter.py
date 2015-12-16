#!/usr/bin/env python
"""
module to recenter the labels based on the manual shifting of tss signal region.

TODO 
    writing fasta seq module 
    how many jobs to packed to a worker node
    score producing module different 
"""

import uuid
import numpy
from Bio import SeqIO

from utils import compressed_pickle


import libpyjobrunner as pg


def reduce_result_score(flat_result):
    """
    collect the results from different worker
    """

    merged_arr = numpy.zeros((len(flat_result), len(flat_result[0])))
    
    for idx, element in enumerate(flat_result):
        merged_arr[idx] = element

    return merged_arr


def reduce_result(flat_result_list):
    """
    collect the results from different worker
    """
    
    merged_arr = numpy.zeros((300, 1000)) ## 3K examples +/-500 center_offset 
    idx = 0 

    for flat_result in flat_result_list:
        for element in flat_result:
            merged_arr[idx] = element
            idx += 1 

    return merged_arr
   

def predict_around_region(args_list):
    """
    predict along the defined region of signal site 
    """

    start_scan, stop_scan, model, data_set = args_list 

    trimed_seq_list = [] 
    for datum in data_set:
        tss_score = numpy.zeros(1000) ## predefined flanking region of 500 nts 

        ## predicting the upstream and downstream of the signal site 
        for idy, shifted_cen in enumerate(xrange(start_scan, stop_scan)):
            model.param["center_pos"] = shifted_cen
            out = model.predict(datum)
            tss_score[idy] = out[0] 

        trimed_seq_list.append(tss_score)
    
    return trimed_seq_list


def load_data_from_fasta(signal, org, data_path, label_type="pos"):
    """
    load examples and labels from fasta file
    """
    
    fn_pos = "%s/%s_sig_%s_example.fa" % (data_path, signal, label_type)
    print "loading: \n %s" % (fn_pos) 

    # load the seq records from fasta file 
    xt_pos = [str(rec.seq) for rec in SeqIO.parse(fn_pos, "fasta")]
    examples =  xt_pos

    print("organism: %s, signal %s,\t num_examples %i" %  (org, signal, len(examples)))

    ret = {"examples": numpy.array(examples)}

    return ret


def load_svm_model(filename):
    """
    load the model from file 
    """
    import bz2 
    import cPickle 

    fh = bz2.BZ2File(filename, "rb")
    model = cPickle.load(fh) 
    fh.close() 

    return model


def signal_center_offset_score(svm_file, org, signal="tss", data_path="SRA-seq"):
    """
    manually look at the position around the original position 
    """
    prefix_fname = "stringtie_dmel_300ex_500bp" 
    center_offset = 500 ## flanking region around signal site 

    ## loading data from fasta file
    data = load_data_from_fasta(signal, org, data_path)
    assert(len(data["examples"]) == len(data["labels"]))

    model_obj = load_svm_model(svm_file)
   
    ## getting the model information 
    center_pos = model_obj.param["center_pos"]
    print "model - center pos: %i, center reg: %i" % (center_pos, center_offset) 

    start_scan = center_pos-center_offset
    stop_scan = center_pos+center_offset

    cnt = 0
    data_set = [] 
    argument_list = [] 

    ## get the individual examples to recenter the signal position manually
    for idx, single_example in enumerate(data["examples"]):  

        datum = [single_example]
        cnt += 1

        if cnt % 50 == 0: ## packing 10 seq to one job 
            data_set.append(datum)

            arg = [start_scan, stop_scan, model_obj, data_set]
            argument_list.append(arg)

            data_set = [] 
            if cnt == 300:
                break
        else:
            data_set.append(datum)
                
    local = False 
    #cluster_resource = {'pvmem':'6gb', 'pmem':'6gb', 'mem':'6gb', 'vmem':'6gb','ppn':'1', 'nodes':'1', 'walltime':'8:00:00'}
    cluster_resource = {'mem':'6000', 'nodes':'1', 'walltime':'08:00'}

    intm_ret = pg.pg_map(predict_around_region, argument_list, param=cluster_resource, local=local, maxNumThreads=1, mem="6000") 
    print "Done with computation"
    
    pred_out_val = reduce_result_score(intm_ret) 
    print "Done reducing the results"

    ## save the scores 
    fname = "%s_%s_score_%s" % (signal, prefix_fname, uuid.uuid1()) 
    compressed_pickle.save(fname, pred_out_val) 
    print( "saving the score in file %s" % fname )


def main():

    org_code  = "H_sapiens"
    model_obj = "tss_04e01036-4ce7-11e5-a179-90b11c0884ee"

    import sys 
    
    try:
        org_code = sys.argv[1]
        model_obj = sys.argv[2]
    except:
        exit(__doc__)

    signal_center_offset_score(model_obj, org_code)
 

if __name__ == "__main__":
    main()

