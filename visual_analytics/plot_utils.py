
import numpy 
import pandas 
import matplotlib.pyplot as plt

import bz2 
import cPickle 

from scipy.signal import argrelextrema


def bar_chart_auroc(file)
    """
    plot a bar chart based on the test performance from diff CV
    """

    file = "hs_s32_trsk_tss.pickle"

    ## loading the data 
    data = data_process(file)

    method = 'individual'
    org_code = 'H_sapiens'
    numb_cvs =  len(data[method][0][org_code])

    ## fetching the best performance based on the validation and test 
    best_eval = numpy.argmax(data[method][0][org_code], axis=1)
    test_perf = numpy.zeros(numb_cvs) 
    for idx, v_idx in enumerate(best_eval):
        test_perf[idx] = data[method][1][org_code][idx][v_idx]

    ## plot settings 
    import pylab
    pylab.figure(figsize=(3,4))
    pylab.rcParams.update({'figure.autolayout':True})

    offset = 0 
    width = 0.10
    seperator = 0.10 

    xlocations = [] 
    min_max = [] 
    rects = [] 

    x_axis_marks = []
    for idk, ele in enumerate(test_perf):
        xlocations.append(offset + 0.05)
        min_max.append(ele) 
        rects.append(pylab.bar(offset, ele, width, color="#88aa33", edgecolor="white"))
        offset += width

        idk += 1 
        cvname = "fold_%d" % idk
        x_axis_marks.append(cvname) 
        
    ## mean bar 
    rects.append(pylab.bar(offset, test_perf.mean(), width, color="#ff9999", edgecolor="white")) 
    xlocations.append(offset + 0.05)
    offset += width

    x_axis_marks.append('mean')

    ## adjusting the axis range 
    min_max.sort()
    ymax = min_max[-1]*1.05 
    ymin = min_max[0]*0.98

    ## set the ticks
    tick_step = 0.02 
    ticks = [tick_step*i for i in xrange(int(round(ymax/tick_step)+1))]

    pylab.yticks(ticks)
    pylab.xticks(xlocations, x_axis_marks, rotation="vertical")

    pylab.xlim(0, offset) 
    pylab.ylim(ymin, ymax) 

    plot_title = "title of the figure" 
    pylab.title(plot_title)

    pylab.gca().get_yaxis().tick_left()
    pylab.gca().get_xaxis().tick_bottom()

    pylab.gca().get_yaxis().grid(True)
    pylab.gca().get_xaxis().grid(False)

    fout_name ="test_performance_with_CV.pdf"
    pylab.savefig(fout_name)


def data_process(fname):
    """
    data loading module
    """

    fh = bz2.BZ2File(fname, "rb")
    data = cPickle.load(fh)
    fh.close()
    
    return data 


def bar_chart():
    """
    draw a bar chart
    """

    ## load data 
    file = "tss_pos_1k_score_70co_ce5b9c8a-4f41-11e5-add3-90e2ba3a73f4"
    tss_score = data_process(file) 
    
    ## max prediction output value from each example
    max_pred_out = numpy.zeros(len(tss_score[0]))
    for idx, ele in enumerate(tss_score):
        pos = numpy.argmax(ele)
        max_pred_out[pos] += 1

    ## position based frequency of examples
    pred_out = numpy.zeros((len(max_pred_out), 2))
    for idx, numb in enumerate(max_pred_out):
        pred_out[idx] = numpy.array([idx, numb])

    df_pred_score = pandas.DataFrame(pred_out)
    
    ## plotting settings  
    fig = plt.figure()
    width = .1
    ind = numpy.arange(len(tss_score[0]))

    plt.bar(ind, df_pred_score[1], color="green", edgecolor='white')
    plt.xlabel('Max prediction output value of positive examples')
    plt.ylabel("frequency")

    fout_name ="max_pred_out_score_examples.pdf"
    plt.savefig(fout_name)


def signal_line(file):
    """
    signal distribution score over positions 
    """

    ## load data 
    file = "tss_score_e7a0db3e-4d4c-11e5-ae67-90b11c0884ee"
    tss_score = data_process(file) 

    ## plot settings
    dt = 1.0
    x_axis_pos = 100
    t = numpy.arange(0, x_axis_pos, dt)

    ## multiple signal lines with different colors 
    plt.plot(t, tss_score[7], 'b-')
    plt.plot(t, tss_score[8], 'g-')
    plt.plot(t, tss_score[9], 'r-')
    plt.plot(t, tss_score[10], 'y-')
    
    ## multiple signal lines with same color 
    #for idx, content in  enumerate(tss_score):
    #    if idx ==10:
    #        break 
    #    plt.plot(t, content, 'y-')
 
    plt.xlim(0.0,100.0)
    plt.xlabel('50 nucleotides upstream and downstream flanking region of tss signal')
    plt.ylabel('model prediction output score')
    plt.grid(True)

    ## annotate the image with localmaxima 
    #plt.annotate('local max', xy=(2, 1), xytext=(3, 1.5),
    #            arrowprops=dict(facecolor='black', shrink=0.05),
    #                        )
    #argrelextrema(tss_score[0], numpy.greater)

    fout_name = "examples_pred_out_score.pdf"
    plt.savefig(fout_name)

