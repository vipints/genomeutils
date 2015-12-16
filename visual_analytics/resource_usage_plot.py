#!/usr/bin/env python 
"""
usage of different resources time mem for different methods

use this script on garage machine to create small legends 
"""

import numpy 
import matplotlib.pyplot as plt 
import matplotlib.path as mpath 


def visualize_memory_usage(out_pdf_file):
    """
    memory usage of each method 
    """

    star = mpath.Path.unit_regular_star(5) 
    circle = mpath.Path.unit_circle()

    verts = numpy.concatenate([star.vertices, circle.vertices[::-1, ...]])
    codes = numpy.concatenate([star.codes, circle.codes]) 

    cut_star = mpath.Path(verts, codes) 

    ## data import FIXME read the csv format file with experiment details  
    ind_mem_taken = numpy.array([715.4726563, 715.6367188, 1413.761719, 1651.101563, 2233, 3333.40625])
    union_mem_taken = numpy.array([718.9804688, 715.5585938, 719.3164063, 2136.8125, 2840.890625, 4212.792969])
    mtl_mem_taken = numpy.array([1232.699219, 1303.710938, 1857.699219, 2365.992188, 3074.230469, 4246.414063])
    mtmkl_mem_taken = numpy.array([2636.988281, 2904.535156, 3258.386719, 3970.226563, 4674.871094, 5847.050781])

    plt.plot(ind_mem_taken, color='#BDB76B', marker=cut_star, markersize=7, linestyle='--')
    plt.plot(union_mem_taken, color='#9370DB', marker=cut_star, markersize=7, linestyle='--')
    plt.plot(mtl_mem_taken, color='#8FBC8B', marker=cut_star, markersize=7, linestyle='--')
    plt.plot(mtmkl_mem_taken, color='#A0522D', marker=cut_star, markersize=7, linestyle='--')

    plt.xticks()

    #TODO the labeling of x axis with # of experiments 
    x_axis = [0, 1, 2, 3, 4, 5]
    labels = [400, 2000, 4000, 8000, 12000, 20000] 
    plt.xticks(x_axis, labels, rotation='vertical') 

    ax = plt.gca() 
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(8)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(8)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    ax.get_xaxis().grid(False)
    ax.get_yaxis().grid(True)

    plt.legend(('individual', 'union', 'mtl', 'mtmkl'), loc='upper left', fontsize=8)
    plt.title("Resource requirement for different methods - memory in megabytes", fontsize=8)
    plt.ylabel('memory in megabyte', fontsize=7)
    plt.xlabel('16 different organisms, each has different set of examples (60% for training the model) used for each method with 1:3 positive, negative combo.', fontsize=7)

    plt.savefig(out_pdf_file)


def time_usage_each_method(out_pdf_file):
    """
    time usage of each method 
    """

    star = mpath.Path.unit_regular_star(5) 
    circle = mpath.Path.unit_circle()

    verts = numpy.concatenate([star.vertices, circle.vertices[::-1, ...]])
    codes = numpy.concatenate([star.codes, circle.codes]) 

    cut_star = mpath.Path(verts, codes) 

    ## FIXME data import part with a csv reader 
    ind_time_taken = numpy.array([56.12435293, 65.80995893, 112.9028471, 179.1727109, 225.9667561, 345.3043289])
    union_time_taken = numpy.array([70.06569886, 100.051934, 103.2509151, 136.4559519, 270.7133741, 304.8211379])
    mtl_time_taken = numpy.array([124.3881071, 355.7790811, 740.989918, 1197.821371, 1529.896782, 2454.422433])
    mtmkl_time_taken = numpy.array([1196.226963, 2115.524649, 3292.124688, 4480.232922, 5862.847364, 8826.01638])

    plt.plot(numpy.log10(ind_time_taken), color='#BDB76B', marker=cut_star, markersize=8, linestyle='--')
    plt.plot(numpy.log10(union_time_taken), color='#9370DB', marker=cut_star, markersize=8, linestyle='--')
    plt.plot(numpy.log10(mtl_time_taken), color='#8FBC8B', marker=cut_star, markersize=8, linestyle='--')
    plt.plot(numpy.log10(mtmkl_time_taken), color='#A0522D', marker=cut_star, markersize=8, linestyle='--')

    plt.xticks()
    x_axis = [0, 1, 2, 3, 4, 5]
    labels = [400, 2000, 4000, 8000, 12000, 20000] 
    plt.xticks(x_axis, labels, rotation='vertical') 
    
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(8)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(8)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    ax.get_xaxis().grid(False)
    ax.get_yaxis().grid(True)

    # FIXME the title and labels for axis to be user provided
    plt.legend(('individual', 'union', 'mtl', 'mtmkl'), loc='upper left', fontsize=8)
    plt.title("Resource requirement for different methods - time in seconds", fontsize=8)
    plt.ylabel('time in seconds (log base 10)', fontsize=7)
    plt.xlabel('16 different organisms, each has different set of examples (60% for training the model) used for each method with 1:3 positive, negative combo.', fontsize=7)

    plt.savefig(out_pdf_file)

 
