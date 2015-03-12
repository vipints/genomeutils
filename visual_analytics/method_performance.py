"""
detailed bar plot to measure the performance of different methods 

Usage: 

from visual_analytics import method_performance as mp 
fname = '29_org_base_0.5_acc.pickle'
organisms, diff_methods, perfomance = mp.data_process(fname) 
mp.detailed_barplot(perfomance, diff_methods, organisms, 'acc.pdf', 'acceptor splice site') 

Requirement:
    pylab: 
    numpy 
"""

from __future__ import division
import numpy 
from collections import defaultdict

def detailed_barplot(data, methods, labels, res_file, plot_title="", ylabel="auROC"):
    """
    visualizing the experiment result of different methods 
    
    @args data: data to plot, for each organism different methods and its performances
    @type data: <defaultdict<org_name:[(method, list)]>
    @args methods: different experiment methods 
    @type methods: list
    @args labels: X axis labels, organism names 
    @type labels: list 
    @args res_file: result file name 
    @type res_file: str
    """
    import pylab 

    #pylab.figure(figsize=(5, 10)) # custom form 
    pylab.figure(figsize=(len(labels), (len(labels)/8)*5)) # 40, 10 # for 10 organisms 
    pylab.rcParams.update({'figure.autolayout': True}) # to fit the figure in canvas 

    width = 0.20
    separator = 0.15
    offset = 0
    num_methods = len(methods)
    
    used_colors = ["#88aa33", "#9999ff", "#ff9999", "#34A4A8"]
    xlocations = []
    
    min_max = [] 
    mean_perf = defaultdict(list) 

    for org_name, details in data.items():
        offset += separator
        rects = [] 

        #xlocations.append(offset + (width*(num_methods*7+2))/2)
        xlocations.append(offset + (width*(num_methods*1))/3)

        for idx, bundles in enumerate(details):
            method, perfs = bundles 
            ##print '\t', method

            best_c = [] 
            for method_perf in perfs: 
                best_c.append(method_perf)

            best_c.sort() 
            min_max.append(best_c[-1])
            mean_perf[method].append(best_c[-1]) # best c over organisms on each method 

            rects.append(pylab.bar(offset, best_c[-1], width, color=used_colors[idx], edgecolor='white'))
            offset += width 

        #offset += separator
        #break
   
    # average of each methods  
    rects_avg = [] 
    offset += separator
    #for meth_color, method_avg in mean_perf.items():
    #    print meth_color
    xlocations.append(offset + (width*(num_methods*1))/3)
    rects_avg.append(pylab.bar(offset, round(sum(mean_perf['union'])/len(labels), 2), width, color = used_colors[0], edgecolor='white'))
    offset += width 
    rects_avg.append(pylab.bar(offset, round(sum(mean_perf['individual'])/len(labels), 2), width, color = used_colors[1], edgecolor='white'))
    offset += width 
    rects_avg.append(pylab.bar(offset, round(sum(mean_perf['mtl'])/len(labels), 2), width, color = used_colors[2], edgecolor='white'))
    offset += width 
    #rects_avg.append(pylab.bar(offset, sum(mean_perf['mtmkl'])/len(labels), width, color = used_colors[3], edgecolor='white'))
    #offset += width 
    
    print 'mean', round(sum(mean_perf['union'])/len(labels), 2), round(sum(mean_perf['individual'])/len(labels), 2), round(sum(mean_perf['mtl'])/len(labels), 2) 

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


def data_process(filename):
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
    # organme - method - perfomance measure 
    data_mat = defaultdict(list) 

    for method, org_perf in myobj.items():
        methods.append(method) 
        
        #for name, perf_meas in org_perf.items():
        for name, perf_meas in org_perf[0].items():

            # organme - method - mean of perfomance measure from different cross validation
            data_mat[name].append((method, perf_meas.mean(axis=0)))
            
    fh.close()
    return data_mat.keys(), methods, data_mat 


def mean_plot_diff_run(data_dir, res_file, signal="cleave signal"):
    """
    visualizing mean performance from different experiments

    @args data_dir: experiment path 
    @type data_dir: str 
    @args alpha: different parameter value 
    @type alpha: list 
    """
    
    import os 
    import re 

    file_mean = defaultdict(list) 
    for file in os.listdir(data_dir):
        #FIXME adjust the next line to get the name from the file name  
        #prefix = re.search('5.org_(\d+)_cleave.pickle', file).group(1) 
        prefix = re.search('5.org_(.+)_cleave.pickle', file).group(1) 
        
        mean_perf = defaultdict(list) 
        organisms, diff_methods, perfomance = data_process("%s/%s" % (data_dir, file)) 

        for org_name, details in perfomance.items():
            for idx, bundles in enumerate(details):
                method, perfs = bundles 
                #print '\t', method

                best_c = [] 
                for method_perf in perfs: 
                    best_c.append(method_perf)

                best_c.sort() 
                #min_max.append(best_c[-1])
                mean_perf[method].append(best_c[-1]) # best c over organisms on each method 
        
        file_mean[prefix].append(sum(mean_perf['union'])/len(organisms)) 
        file_mean[prefix].append(sum(mean_perf['individual'])/len(organisms)) 
        file_mean[prefix].append(sum(mean_perf['mtl'])/len(organisms)) 

    ## plotting function 
    import pylab 

    pylab.figure(figsize=(5, 10)) # custom form 
    #pylab.figure(figsize=(len(labels), (len(labels)/8)*5)) # 40, 10 # for 10 organisms 
    pylab.rcParams.update({'figure.autolayout': True}) # to fit the figure in canvas 

    width = 0.20
    separator = 0.15

    offset = 0
    num_methods = len(diff_methods)
    
    used_colors = ["#88aa33", "#9999ff", "#ff9999", "#34A4A8"]
    xlocations = []
    min_max = [] 
    for parameter, mean_meth_perf in sorted(file_mean.items()):
        offset += separator
        rects = [] 

        xlocations.append(offset + (width*(num_methods*1))/3)
        for idx, nb in enumerate(mean_meth_perf):
            min_max.append(nb)
            rects.append(pylab.bar(offset, nb, width, color=used_colors[idx], edgecolor='white'))
            offset += width 

    offset += separator
            
    min_max.sort() 
    ymax = min_max[-1]*1.1
    ymin = min_max[0]*0.9

    tick_step = 0.05
    ticks = [tick_step*i for i in xrange(round(ymax/tick_step)+1)]
    pylab.yticks(ticks)
    labels = ['df=%s' % x for x in sorted(file_mean.keys())]
    pylab.xticks(xlocations, labels, rotation="vertical") 

    fontsize=17
    ax = pylab.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    pylab.xlim(0, offset)
    pylab.ylim(ymin, ymax)
    
    plot_title = signal 
    pylab.title(plot_title)
    pylab.gca().get_xaxis().tick_bottom()
    pylab.gca().get_yaxis().tick_left()

    pylab.gca().get_yaxis().grid(True)
    pylab.gca().get_xaxis().grid(False)
    pylab.legend(tuple(rects), tuple(diff_methods))

    ylabel = "auROC"
    pylab.ylabel(ylabel, fontsize = 15)
    pylab.savefig(res_file) 

