"""
Usage: 

Requirement:
    pylab: 
"""

from __future__ import division
from collections import defaultdict

def barplot(eval_data, test_data, methods, labels, res_file, plot_title="", ylabel="auROC"):
    """
    visualizing the experiment result of different methods 
    
    @args eval_data: data to plot, for each organism different methods and its performances
    @type eval_data: <defaultdict<org_name:[(method, list)]>
    @args test_data: data to plot, for each organism different methods and its performances
    @type test_data: <defaultdict<org_name:[(method, list)]>
    @args methods: different experiment methods 
    @type methods: list
    @args labels: X axis labels, organism names 
    @type labels: list 
    @args res_file: result file name 
    @type res_file: str
    """

    import numpy 
    import pylab 

    pylab.figure(figsize=(10, 10)) # custom form 
    pylab.rcParams.update({'figure.autolayout': True}) # to fit the figure in canvas 

    width = 0.20
    separator = 0.15
    offset = 0
    num_methods = len(methods)
    used_colors = ["#88aa33", "#9999ff", "#ff9999", "#34A4A8"]
    xlocations = []
    
    min_max = [] 
    mean_perf = defaultdict(list) 
    index_method = dict( individual = 0,
                        union = 1,
                        mtl = 2, 
                        mtmkl = 3)

    for org_name, details in eval_data.items():
        offset += separator
        rects = [] 
        xlocations.append(offset + (width*(num_methods*1))/3)

        for idx, bundles in enumerate(details):
            method, perfs = bundles 
            print '\t', method
            best_c = numpy.argmax(perfs)
            best_c_score = test_data[org_name][index_method[method]][1][best_c] 

            min_max.append(best_c_score)
            mean_perf[method].append(best_c_score)

            rects.append(pylab.bar(offset, best_c_score, width, color=used_colors[idx], edgecolor='white'))
            offset += width 

    # average of each methods  
    rects_avg = [] 
    offset += separator
    xlocations.append(offset + (width*(num_methods*1))/3)

    rects_avg.append(pylab.bar(offset, sum(mean_perf['individual'])/len(labels), width, color = used_colors[0], edgecolor='white'))
    offset += width 

    rects_avg.append(pylab.bar(offset, sum(mean_perf['union'])/len(labels), width, color = used_colors[1], edgecolor='white'))
    offset += width 

    rects_avg.append(pylab.bar(offset, sum(mean_perf['mtl'])/len(labels), width, color = used_colors[2], edgecolor='white'))
    offset += width 
    
    rects_avg.append(pylab.bar(offset, sum(mean_perf['mtmkl'])/len(labels), width, color = used_colors[3], edgecolor='white'))
    offset += width 

    print 'individual - mean perf', round(sum(mean_perf['individual'])/len(labels), 2)
    print 'union - mean perf', round(sum(mean_perf['union'])/len(labels), 2)
    print 'mtl - mean perf', round(sum(mean_perf['mtl'])/len(labels), 2)
    print 'mtmkl - mean perf', round(sum(mean_perf['mtmkl'])/len(labels), 2) 

    offset += separator
    labels.append('Mean')

    # determine the extreams 
    min_max.sort() 
    ymax = min_max[-1]*1.1
    ymin = min_max[0]*0.9

    # set ticks
    tick_step = 0.05
    ticks = [tick_step*i for i in xrange(round(ymax/tick_step)+1)]

    pylab.yticks(ticks)
    pylab.xticks(xlocations, labels, rotation="vertical") 

    fontsize=17
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    pylab.xlim(0, offset)
    pylab.ylim(ymin, ymax)
    
    pylab.title(plot_title)
    pylab.gca().get_xaxis().tick_bottom()
    pylab.gca().get_yaxis().tick_left()

    pylab.gca().get_yaxis().grid(True)
    pylab.gca().get_xaxis().grid(False)

    #pylab.legend(tuple(rects[0:21:7]), tuple(methods))
    pylab.legend(tuple(rects), tuple(methods))

    pylab.ylabel(ylabel, fontsize = 15)
    pylab.savefig(res_file) 


def parse_data(filename):
    """
    process the pickle file to get data frames

    @args filename: pickle file from experiment run 
    @type filename: bz2 pickle file 
    """
    
    import bz2 
    import cPickle 

    fh = bz2.BZ2File(filename, 'rb')
    myobj = cPickle.load(fh) 

    methods = [] 
    eval_perf_tmp = defaultdict(list) 
    test_perf_tmp = defaultdict(list) 

    for method, org_perf in myobj.items():
        methods.append(method) 
        
        for name, perf_meas in org_perf[0].items():## eval performance  
            eval_perf_tmp[name].append((method, perf_meas.mean(axis=0)))# organme - method - mean of perfomance measure from different cross validation
        for name, test_meas in org_perf[1].items():## test performance 
            test_perf_tmp[name].append((method, test_meas.mean(axis=0)))
            
    fh.close()

    diff_methods = ['individual', 'union', 'mtl', 'mtmkl'] ## pre-defined methods for learning techniques 
    assert (set(methods)==set(diff_methods)), "methods from pickle file %s != %s" % (methods, diff_methods)

    # making an order for the experiments
    eval_perf = defaultdict(list)
    test_perf = defaultdict(list)
    for order_meth in diff_methods:
        for org in eval_perf_tmp.keys():
            for methods in eval_perf_tmp[org]:
                if methods[0] == order_meth:
                    eval_perf[org].append(methods)

    for order_meth in diff_methods:
        for org in test_perf_tmp.keys():
            for methods in test_perf_tmp[org]:
                if methods[0] == order_meth:
                    test_perf[org].append(methods)

    return eval_perf.keys(), diff_methods, eval_perf, test_perf 

