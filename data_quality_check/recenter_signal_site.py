#!/usr/bin/env python
"""
module to recenter the genomic signal sequence based on the manual shifting.
"""
## standard modules 
import uuid
import numpy
import random 
from Bio import SeqIO

## custom modules
from utils import compressed_pickle
from promoter_kernel import ShogunPredictor

## distributed computing 
import libpyjobrunner as pg


def coshuffle(*args):
    """
    will shuffle target_list and apply
    same permutation to other lists

    >>> coshuffle([2, 1, 3], [4, 2, 8], [6, 3, 12])
    ([5, 3, 2, 1, 4], [5, 3, 2, 1, 4], [5, 3, 2, 1, 4])
    """ 

    assert len(args) > 0, "need at least one list"
    num_elements = len(args[0])

    for arg in args:
        assert len(arg) == num_elements, "length mismatch"

    idx = range(num_elements)
    random.shuffle(idx)
    new_lists = []

    for arg in args:
        new_lists.append([arg[i] for i in idx])

    return tuple(new_lists)


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

    print("organism: %s, signal %s,\t num_labels: %i,\t num_examples: %i" %  (org, signal, len(labels), len(examples)))

    examples_shuffled, labels_shuffled = coshuffle(examples, labels)
    ret = {"examples": numpy.array(examples_shuffled), "labels": numpy.array(labels_shuffled)}

    return ret


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


def predict_site_region(args_list):
    """
    nucleotide level prediction around the signal center region
    """

    start_scan, stop_scan, model, data_set = args_list 

    ex_seq_list = [] 
    cen_flank_region = stop_scan - start_scan

    for datum in data_set:
        tss_score = numpy.zeros(cen_flank_region) ## center flanking region 

        ## predicting the near regions of the signal, 
        for idy, shifted_cen in enumerate(xrange(start_scan, stop_scan)):
            model.param["center_pos"] = shifted_cen
            out = model.predict(datum)
            tss_score[idy] = out[0] 

        ex_seq_list.append(tss_score)
    
    return ex_seq_list


def reduce_pred_score(flat_result_list):
    """
    collect the scores from different worker
    """
    
    idx = 0 
    merged_arr = numpy.zeros((len(flat_result_list), len(flat_result_list[0])))

    for flat_result in flat_result_list:
        for element in flat_result:
            merged_arr[idx] = element
            idx += 1 

    return merged_arr


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


def reduce_modified_seq(flat_result):
    """ return the result 
    """
    
    return flat_result
    
    
def write_fasta_rec(seq_list_total, signal, ex_type):
    """ write fasta file from the seq records 
    """

    seq_type = {"pos":"+1", "neg":"-1"}
    fname = "%s_sig_%s_ex.fa" % (signal, ex_type) ## trimmed sequence 

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    out_fh = open(fname, "w") ## write examples in fasta format 
    for seq_list in seq_list_total: 
        for seq_rec in seq_list:
            faseq_out = SeqRecord(Seq(seq_rec), id="%s_ex_%s" % (signal, uuid.uuid1()), description=seq_type[ex_type])
            out_fh.write(faseq_out.format('fasta'))
    out_fh.close() 


def data_process_depot(svm_file, org, example_type, signal, data_path, num_seqs, center_offset):
    """
    get the input data to do computation
    """

    ## loading data
    data = load_examples_from_fasta(signal, example_type, org, data_path)
    assert(len(data["examples"]) == len(data["labels"]))

    ## unpack the model
    model = load_model(svm_file)
   
    ## getting the model information 
    center_pos = model.param["center_pos"]
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

        if cnt % num_seqs == 0: ## packing 10 seq to one job 
            data_set.append(datum)

            arg = [start_scan, stop_scan, model, data_set]
            argument_list.append(arg)

            data_set = [] 
            if cnt == 2:
                break 
        else:
            data_set.append(datum)
    
    return argument_list


def shift_signal_position(svm_file, org, example_type="pos", signal="tss", data_path="SRA-rnaseq"):
    """
    manually look at the position around the original position 
    """

    local = False ## switch between local and compute cluster 
    ## cluster compute options   
    cluster_resource = {'pvmem':'5gb', 'pmem':'5gb', 'mem':'5gb', 'vmem':'5gb','ppn':'1', 'nodes':'1', 'walltime':'24:00:00'}

    num_seq_ex = 2 ## number of sequences are in a single job  
    center_offset = 50 ## nearby regions 
    args_req_list = data_process_depot(svm_file, org, example_type, signal, data_path, num_seq_ex, center_offset)

    ## job dispatching 
    intm_ret = pg.pg_map(recenter_examples, args_req_list, param=cluster_resource, local=local, maxNumThreads=1, mem="5gb") 
    print("Done with trimming example sequences")

    fixed_example_seq = reduce_modified_seq(intm_ret) 
    print("Done with collecting the trimmed examples")
        
    write_fasta_rec(fixed_example_seq, signal, example_type) 
    print("Done with writing examples in fasta format")


def calculate_pred_score(svm_file, org, example_type="pos", signal="tss", data_path="SRA-rnaseq"):
    """
    calculate svm prediction score around the true signal site
    """

    local = False ## switch between local and compute cluster 
    ## cluster compute options   
    cluster_resource = {'pvmem':'8gb', 'pmem':'8gb', 'mem':'8gb', 'vmem':'8gb','ppn':'1', 'nodes':'1', 'walltime':'24:00:00'}
    #cluster_resource = {'mem':'6000', 'nodes':'1', 'walltime':'08:00'}

    num_seq_ex = 10 ## number of sequences are in single job 
    center_offset = 500 ## nearby regions  
    args_req_list = data_process_depot(svm_file, org, example_type, signal, data_path, num_seq_ex)

    intm_ret = pg.pg_map(predict_site_region, args_req_list, param=cluster_resource, local=local, maxNumThreads=1, mem="8gb") 
    print("Done with calculating the score for center region of example sequences")

    pred_out_val = reduce_pred_score(intm_ret) 
    print("Done with collecting scores from different workers")

    ## save the scores 
    fname = "%s_%s_ex_pred_score_%s" % (signal, example_type, uuid.uuid1()) 
    compressed_pickle.save(fname, pred_out_val) 

    print("saving the scores in file %s" % fname)


def main():

    import sys 

    try:
        org_code = sys.argv[1]
        svm_file_name = sys.argv[2]
        data_location = sys.argv[3]
    except:
        exit(__doc__)

    signal = "tss"
    example_type = "pos"

    shift_signal_position(svm_file_name, org_code, example_type, signal, data_location)


if __name__ == "__main__":
    main()
