from __future__ import division
import numpy 
import matplotlib.pyplot as plt

def plot_corrcf(rho_mat, x_label, y_label, out_file):
    """
    plotting the correlation matrix

    @args rho_mat: a matrix with correlation coefficients 
    @type rho_mat: 2-D numpy array 
    @args x_label, y_label: X-Y axis label names 
    @type x_label, y_label: 1-D numpy array 
    @args out_file: name of destination file
    @type out_file: str 
    """
    
    # adjusting the canvas size
    fig = plt.figure(figsize=(len(x_label)/2.667, len(y_label)/2.667))

    plt.pcolor(rho_mat)
    plt.colorbar()
    
    # adjusting the axis labels 
    plt.xticks(numpy.arange(0.5, len(x_label)+0.5), x_label, rotation='vertical')
    plt.yticks(numpy.arange(0.5, len(y_label)+0.5), y_label) 

    # save the figure
    plt.savefig(out_file, bbox_inches='tight')

