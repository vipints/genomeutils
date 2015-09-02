
import numpy as np
import matplotlib.pyplot as plt

import bz2 
import cPickle 

from scipy.signal import argrelextrema



def data_process(fname):
    """
    data loading module
    """

    fh = bz2.BZ2File(fname, "rb")
    data = cPickle.load(fh)
    fh.close()
    
    return data 



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
    t = np.arange(0, x_axis_pos, dt)

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
    #argrelextrema(tss_score[0], np.greater)

    fout_name = "examples_pred_out_score.pdf"
    plt.savefig(fout_name)

