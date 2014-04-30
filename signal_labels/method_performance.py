"""
detailed bar plot to measure the performance of different methods 
"""
from __future__ import division
from collections import defaultdict
import numpy 

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

    #pylab.figure(figsize=(len(labels)*4, (len(labels)/10)*5)) # 40, 10 # for 10 organisms 
    pylab.figure(figsize=(len(labels), (len(labels)/8)*5)) # 40, 10 # for 10 organisms 
    pylab.rcParams.update({'figure.autolayout': True}) # to fit the figure in canvas 

    width = 0.20
    separator = 0.15

    offset = 0
    num_methods = len(methods)
    
    used_colors = ["#88aa33", "#9999ff", "#ff9999"]
    xlocations = []
    
    min_max = [] 
    mean_perf = defaultdict(list) 

    for org_name, details in data.items():
        
        offset += separator
        rects = [] 
        
        print 'organism', org_name

        #xlocations.append(offset + (width*(num_methods*7+2))/2)
        xlocations.append(offset + (width*(num_methods*1))/3)

        for idx, bundles in enumerate(details):
            method, perfs = bundles 

            best_c = [] 
            for method_perf in perfs: 
                best_c.append(method_perf)

            best_c.sort() 
            min_max.append(best_c[-1])
            mean_perf[(method, used_colors[idx])].append(best_c[-1]) # average of best c over organisms on each method 

            rects.append(pylab.bar(offset, best_c[-1], width, color=used_colors[idx], edgecolor='white'))
            offset += width 

        offset += separator
        #break
   
    # average of each methods  
    rects_avg = [] 
    offset += separator
    for meth_color, method_avg in mean_perf.items():
        rects_avg.append(pylab.bar(offset, sum(method_avg)/len(labels), width, color = meth_color[1], edgecolor='white'))
        offset += width 
    offset += separator

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
        
        for name, perf_meas in org_perf.items():

            # organme - method - mean of perfomance measure from different cross validation
            data_mat[name].append((method, perf_meas.mean(axis=0)))
            
    fh.close()
    return data_mat.keys(), methods, data_mat 

"""
from signal_labels import method_performance
fname = '29_org_base_0.5_acc.pickle'
organisms, diff_methods, perfomance = method_performance.data_process(fname) 
method_performance.detailed_barplot(perfomance, diff_methods, organisms, 'acc.pdf', 'acceptor splice site') 
"""