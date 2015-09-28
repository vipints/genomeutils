#!/usr/bin/env python 
"""
setting up experiments for tss signal using string kernels with promoter features
"""

import time
import numpy as np
from itertools import chain
import uuid
from sklearn.grid_search import ParameterGrid

from multitask_cross_validation import MultitaskShuffleSplitThreeWay, subset_mtl_data
from model_mkl import ModelIndividual

from data_processing import data_loader, helper 

import libpyjobrunner as pg


def reduce_result(list_kv_pairs):
    """
    list of key-value pairs
    """

    task_keys, fold_indices, grid_indices = zip(*[item.keys()[0] for item in list_kv_pairs])
    num_folds = len(set(fold_indices))
    num_grid_points = len(set(grid_indices))
    
    task_keys = list(set(task_keys))

    perf_dev = dict((org, np.zeros((num_folds, num_grid_points))) for org in task_keys)
    perf_test = dict((org, np.zeros((num_folds, num_grid_points))) for org in task_keys)

    # fill return matrix
    for item in list_kv_pairs:
        (task_key, fold_idx, grid_idx) = item.keys()[0]
        perf_dev[task_key][fold_idx, grid_idx] = item.values()[0]["prc_dev"]
        perf_test[task_key][fold_idx, grid_idx] = item.values()[0]["prc_test"]

    return perf_dev, perf_test


def compute_core(arg_list):
    """
    call to be distributed onto the cluster
    """

    signal, method, fold_idx, train_idx, dev_idx, test_idx, grid_idx, grid_point = arg_list

    data = data_loader.load_all(signal)

    data_train = subset_mtl_data(data, train_idx)
    data_dev = subset_mtl_data(data, dev_idx)
    data_test = subset_mtl_data(data, test_idx)

    model = method(cost=grid_point["cost"])
    model.train(data_train)
    
    out_test = model.predict(data_test)
    out_dev = model.predict(data_dev)

    ret = []
    # return a list of key-value pairs
    for key in data_test.keys():
        perf_prc_dev = helper.calcroc(out_dev[key], data_dev[key]["labels"])[0]
        perf_prc_test = helper.calcroc(out_test[key], data_test[key]["labels"])[0]
        reduce_key = (key, fold_idx, grid_idx)
        ret.append({reduce_key : {"prc_dev": perf_prc_dev, "prc_test": perf_prc_test}})

    return ret


def setup_splits(signal, method_name, method, param, num_folds, test_size, random_state):
    """
    splitting the example data into train/test/validation group 
    """

    data = data_loader.load_all(signal)
    sizes = dict((org, len(data[org]["labels"])) for org in data.keys())

    # set up splitting strategy
    kf = MultitaskShuffleSplitThreeWay(sizes, n_iter=num_folds, indices=True, test_size=test_size*2, random_state=random_state)

    param_grid = list(ParameterGrid(param))
    argument_list = []

    for fold_idx, (train_idx, dev_idx, test_idx) in enumerate(kf):
        for grid_idx, grid_point in enumerate(param_grid):
            arg = [signal, method, fold_idx, train_idx, dev_idx, test_idx, grid_idx, grid_point]
            argument_list.append(arg)

    local = False 
    max_num_threads = 2

    if method_name in ['union', 'individual']:
        param = {'vmem':'4gb', 'pvmem':'4gb', 'pmem':'4gb', 'mem':'4gb', 'ppn':'1', 'nodes':'1', 'walltime':'2:00:00'}
        intermediate_ret = pg.pg_map(compute_core, argument_list, param=param, local=local, maxNumThreads=1, mem="4gb")

    #import ipdb 
    #ipdb.set_trace()

    print "DONE with computation"

    flat_intermediate = list(chain.from_iterable(intermediate_ret))
    perf_dev, perf_test = reduce_result(flat_intermediate)

    print "DONE reducing"

    return perf_dev, perf_test



def run_single_org_signal(signal):
    """
    wrapper for experiment for particular signal
    """

    exp_prefix = "arts_set_1_shift32"
    results = {}

    # set up hyper-parameters C
    costs = np.logspace(-2, 2, num=7, endpoint=True) ## FIXME change me back  
    param = {'cost': costs}

    methods = {"individual": ModelIndividual}

    times = {}
    for method_name, method in methods.items():
        t0 = time.time()
        num_folds = 10
        perf_orgs = setup_splits(signal, method_name, method, param, num_folds, 0.2, 42)

        fn_pickle = "%s_%s_%s.pickle" % (exp_prefix, method_name, signal)
        print("saving", fn_pickle)

        helper.save(fn_pickle, perf_orgs)
        results[method_name] = perf_orgs
 
        time_taken = time.time() - t0
        times[method_name] = time_taken
        print "time taken for method", method_name, time_taken

    fn_res = "%s_%s.pickle" % (exp_prefix, signal)
    helper.save(fn_res, results)
    print times 


def main():
    """
    defining experiments for different genomic signal 
    """
    
    signal = "tss" 

    run_single_org_signal(signal)


if __name__=="__main__":
    main()
