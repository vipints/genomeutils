"""
detailed bar plot to measure the performance of different methods 
"""
from __future__ import division
from collections import defaultdict
import numpy 

def detailed_barplot(data, methods, plot_title="", ylabel="auROC"):
    """

    """
    
    import pylab 
    pylab.figure(1,(10,30))

    mod_methods = []
    for l in methods:
        mod_methods.append(l)
    
    width = 0.20
    separator = 0.15

    offset = 0
    num_methods = len(methods)
    
    used_colors = ["green", "blue", "red"]
    xlocations = []

    min_max = [] 
    labels = [] 

    for org_name, details in data.items():
        
        offset += separator
        rects = [] 
        
        print 'organism', org_name

        for idx, bundles in enumerate(details):
            method, perfs = bundles 
            
            xlocations.append(offset + num_methods*width/2)
            labels.append('%s_%s' % (org_name, method))

            for method_perf in perfs: 
                
                min_max.append(method_perf) 
                rects.append(pylab.bar(offset, method_perf, width, color=used_colors[idx]))
                offset += width 
        
            offset += separator
        offset += separator
        break
     
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

    #pylab.legend(tuple(rects), tuple(mod_methods))

    pylab.ylabel(ylabel, fontsize = 15)
    pylab.savefig('f1-test.pdf') 


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
        
        for name, perf_meas in org_perf.items():

            # organme - method - mean of perfomance measure from different cross validation
            data_mat[name].append((method, perf_meas.mean(axis=0)))
            
    fh.close()
    return data_mat.keys(), methods, data_mat 

