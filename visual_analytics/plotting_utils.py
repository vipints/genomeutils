
import numpy 
import pandas 
import matplotlib.pyplot as plt

import bz2 
import cPickle 

from scipy.signal import argrelextrema


def plot_distance_matrix(distance_csv_file, outfile="similarity_matrix_plot.pdf"):
    """
    plotting distance matrix between organisms 

    the distance between the organisms are calculated based on the difference in their sequence composition 
    """

    import pylab 

    distance = pandas.read_csv(distance_csv_file, header=0)
    C = numpy.tril(distance)
    sim = 1-distance
    C = numpy.tril(sim)
    N = sim.shape[1]
    C = numpy.ma.masked_array(C, C == 0)

    A = numpy.array([(y, x) for x in range(N, -1, -1) for y in range(N + 1)])
    t = numpy.array([[0.5, 1], [0.5, -1]])
    A = numpy.dot(A, t)
    X = A[:, 1].reshape(N + 1, N + 1)
    Y = A[:, 0].reshape(N + 1, N + 1)
    fig = pylab.figure(figsize=(20,20))
    ax = fig.add_subplot(121, frame_on=False, aspect=2.0)
    ax.set_xticks([])
    ax.set_yticks([])
    caxes = pylab.pcolormesh(X, Y, np.flipud(C), axes=ax)
    ax.set_xlim(right=0)
    fig.savefig(outfile, bbox_inches='tight')


def bar_chart_auroc(file, fout_name, plot_title="performance"):
    """
    plot a bar chart based on the test performance from diff CV
    """

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
    pylab.figure(figsize=(5,4))
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
        rects.append(pylab.bar(offset, ele, width, color="#73C6B6", edgecolor="white"))
        offset += width

        idk += 1 
        cvname = "fold_%d" % idk
        x_axis_marks.append(cvname) 
        
    ## mean bar 
    rects.append(pylab.bar(offset, test_perf.mean(), width, color="#EC7063", edgecolor="white")) 
    xlocations.append(offset + 0.05)
    offset += width

    x_axis_marks.append('mean')

    ## adjusting the axis range 
    min_max.sort()

    ymax = min_max[-1]*1.03 
    ymin = min_max[0]*0.98

    ## set the ticks
    tick_step = 0.01 
    ticks = [tick_step*i for i in xrange(int(round(ymax/tick_step)+1))]

    pylab.yticks(ticks, fontsize=8)
    pylab.xticks(xlocations, x_axis_marks, rotation="vertical", fontsize=8)

    pylab.xlim(0, offset) 
    pylab.ylim(ymin, ymax) 

    pylab.title(plot_title, fontsize=8)

    pylab.gca().get_yaxis().tick_left()
    pylab.gca().get_xaxis().tick_bottom()

    pylab.gca().get_yaxis().grid(True)
    pylab.gca().get_xaxis().grid(False)

    pylab.ylabel("auROC", fontsize=8)

    pylab.savefig(fout_name)


def data_process(fname):
    """
    data loading module
    """

    fh = bz2.BZ2File(fname, "rb")
    data = cPickle.load(fh)
    fh.close()
    
    return data 


def bar_chart(score_file, fout_name, plot_title="position relative to true TSS"):
    """
    draw a bar chart based on the max. prediction score of each example 
    """

    ## load data 
    tss_score = data_process(score_file) 
    x_axis_row = len(tss_score[0])
    
    ## max prediction output value from each example
    max_pred_out = numpy.zeros(x_axis_row)
    for idx, ele in enumerate(tss_score):
        pos = numpy.argmax(ele)
        max_pred_out[pos] += 1

    ## position based frequency of examples
    pred_out = numpy.zeros((len(max_pred_out), 2))
    for idx, numb in enumerate(max_pred_out):
        pred_out[idx] = numpy.array([idx, numb])

    #col1-index and col2-score 
    df_pred_score = pandas.DataFrame(pred_out)

    ## plotting settings  
    width = 0.1
    fig = plt.figure()

    ind = numpy.arange(x_axis_row)
    #plt.bar(ind, df_pred_score[1], color="#7FB3D5", edgecolor='#7FB3D5')
    plt.bar(ind, df_pred_score[1], color="#73C6B6", edgecolor='#73C6B6')

    tick_step = 100 ## step size defined by the number of xaxis ticks  
    xlocations = [tick_step*i+tick_step for i in xrange(x_axis_row/tick_step)] 

    plt.yticks(fontsize=9)
    plt.xticks(xlocations, fontsize=9) 

    plt.title(plot_title, fontsize=9)
    plt.xlabel('svm max. prediction output value of examples', fontsize=9)
    plt.ylabel("features frequency", fontsize=9)

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
    plt.xlabel('50 nucleotides upstream and downstream flanking region of tss signal', fontsize=4)
    plt.ylabel('model prediction output score')
    plt.grid(True)

    ## annotate the image with localmaxima 
    #plt.annotate('local max', xy=(2, 1), xytext=(3, 1.5),
    #            arrowprops=dict(facecolor='black', shrink=0.05),
    #                        )
    #argrelextrema(tss_score[0], numpy.greater)

    fout_name = "examples_pred_out_score.pdf"
    plt.savefig(fout_name)


def max_pred_out_range(file):
    """
    calculate the range of max pred output from each position on center offset
    """
    from collections import defaultdict

    file = "tss_pos_1k_score_66302e2c-4e66-11e5-9783-90b11c42b574"
    data = data_process(file) 

    max_pred_out = numpy.zeros(len(data[0]))

    max_pos_score_range = defaultdict(list) 

    for idx, ele in enumerate(data):
        pos = numpy.argmax(ele)
        max_pred_out[pos] += 1
        max_pos_score_range[pos].append(ele.max())
        
    pred_max_range = numpy.zeros((len(data[0]), 3))
    for idx, freq in enumerate(max_pred_out):
        max_pos_score_range[idx].sort()
        min_value = max_pos_score_range[idx][0] 
        max_value = max_pos_score_range[idx][-1] 

        pred_max_range[idx] = numpy.array([freq, min_value, max_value])
        
    df_pred_max_range = pandas.DataFrame(pred_max_range)  

    ## bar chart to show frequency in each positions
    ## plotting settings  
    fig = plt.figure()
    width = .1
    ind = numpy.arange(len(data[0]))

    plt.bar(ind, df_pred_max_range[0], color="green", edgecolor='white')
    plt.xlabel('Max prediction output value of positive examples', fontsize=10)
    plt.ylabel("frequency", fontsize=10)

    plt.yticks()

    fout_name ="max_pred_out_score_pos_examples.pdf"
    plt.savefig(fout_name)


