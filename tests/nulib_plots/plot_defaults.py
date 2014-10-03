#/usr/bin/env python

import scipy as sp
import scipy.optimize as spo
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import itertools
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='palatino')
rc('font', weight='bold')
rc('mathtext', default='sf')
rc("lines", markeredgewidth=2)
rc("lines", linewidth=3)
rc('axes', labelsize=30) #24
rc("axes", linewidth=2) #2)
rc('xtick', labelsize=28)
rc('ytick', labelsize=28)
rc('legend', fontsize=20) #16
rc('xtick.major', pad=8) #8)
rc('ytick.major', pad=8) #8)
rc('xtick.major', size=13) 
rc('ytick.major', size=13) 
rc('xtick.minor', size=7) 
rc('ytick.minor', size=7) 

# some colors
green1='#00ff00'
orange='#ffba00'

def setup_plot_with_vertical_panels(fig, left, bottom, heights, width, panels):
    """left: left starting point of panels
       bottom: bottom starting point of panels
       width: width of panels
       hieghts: an array of individual panel heights
    """
    
    # setup plotting rectangles
    rects = []
    total_height = 0
    for ii in range(0,panels):
        current_left = left
        current_bottom = bottom + total_height
        current_height = heights[ii]
        current_width = width
        total_height += current_height
        rects.append([current_left, current_bottom, current_width, current_height])
    
    # Add axes and
    # share x-axis of bottom panel and hide x-axis labeling of upper panels
    ax = []
    ax.append(fig.add_axes(rects[0]))
    for ii in range(1, panels):
        ax.append(fig.add_axes(rects[ii], sharex=ax[0]))
        for label in ax[ii].get_xticklabels():
            label.set_visible(False)

    return rects, ax


def set_tick_sizes(ax, major, minor):
    for l in ax.get_xticklines() + ax.get_yticklines():
        l.set_markersize(major)
    for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
        tick.tick1line.set_markersize(minor)
        tick.tick2line.set_markersize(minor)
    ax.xaxis.LABELPAD=10.
    ax.xaxis.OFFSETTEXTPAD=10.

def set_ticklines(ax,major,minor):
    ticklines = ax.xaxis.get_majorticklines()
    plt.setp(ticklines,mew=major)
    ticklines = ax.xaxis.get_minorticklines()
    plt.setp(ticklines,mew=minor)
    ticklines = ax.yaxis.get_majorticklines()
    plt.setp(ticklines,mew=major)
    ticklines = ax.yaxis.get_minorticklines()
    plt.setp(ticklines,mew=minor)


###################################################################
# various possible tick label formats

def supermongolike(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = s
    else:
        s = '%1.4e ' % value
        snew = ScientificNotation(s)
        snew = "$\mathbf{"+snew+"}$"

    return snew


def ScientificNotation(s):
    # transform 1e+004 into 1e4, for example
    tup = s.split('e')
    mantissa = tup[0].rstrip('0').rstrip('.')
    sign = tup[1][0].replace('+', '')
    exponent = tup[1][1:].lstrip('0')
    res = '%se{%s%s}' %(mantissa, sign, exponent)
    return res.replace('e{}','').replace('e',r'{\times}10^')
        



def supermongolike2(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "$\mathbf{"+s+"}$"
    else:
        s = '%1.4e ' % value
        snew = ScientificNotation2(s)
        snew = "$\mathbf{"+snew+"}$"

    return snew

def supermongolike2a(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "$\mathbf{"+s+"}$"
    else:
        if(value == 1.0e0):
            s = '1'
            snew = "$\mathbf{"+s+"}$"
        else:
            if(value == 10.0e0):
                s = '10'
                snew = "$\mathbf{"+s+"}$"
            else:
                s = '%1.4e ' % value
                snew = ScientificNotation2(s)
                snew = "$\mathbf{"+snew+"}$"
    return snew


def supermongolike3(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "$\mathbf{"+s+"}$"
    else:
        s = '%1.4e ' % value
        snew = s#stripe10(s)
        snew = "$\mathbf{"+snew+"}$"

    return snew

def supermongolike4(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "$\mathbf{"+s+"}$"
    else:
        s = '%2.2f ' % value
        snew = s#stripe10(s)
        snew = "$\mathbf{"+snew+"}$"

    return snew


def supermongolike5(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "$\mathbf{"+s+"}$"
    else:
        s = '%2d ' % value
        snew = s#stripe10(s)
        snew = "$\mathbf{"+snew+"}$"

    return snew

def supermongolike6(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "$\mathbf{"+s+"}$"
    else:
        s = '%2.1f ' % value
        snew = s#stripe10(s)
        snew = "$\mathbf{"+snew+"}$"

    return snew

def supermongolike7(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "$\mathbf{"+s+"}$"
    else:
        s = '%2.2f ' % value
        snew = s#stripe10(s)
        snew = "$\mathbf{"+snew+"}$"

    return snew

def supermongolikeN2(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "${"+s+"}$"
    else:
        if(value == 1.0e0):
            s = '1'
            snew = "${\ 10^0}$"
        else:
            s = '%1.4e ' % value
            snew = ScientificNotation2(s)
            snew = "${"+snew+"}$"

    return snew

def supermongolikeN2a(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "${"+s+"}$"
    else:
        if(value == 1.0e0):
            s = '1'
            snew = "${"+s+"}$"
        else:
            if(value == 10.0e0):
                s = '10'
                snew = "${"+s+"}$"
            else:
                s = '%1.4e ' % value
                snew = ScientificNotation2(s)
                snew = "${"+snew+"}$"
    return snew


def supermongolikeN3(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "${"+s+"}$"
    else:
        s = '%1.4e ' % value
        snew = s#stripe10(s)
        snew = "${"+snew+"}$"

    return snew

def supermongolikeN4(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "${"+s+"}$"
    else:
        s = '%2.2f ' % value
        snew = s#stripe10(s)
        snew = "${"+snew+"}$"

    return snew

def supermongolikeN4a(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "${"+s+"}$"
        return snew
    if(value < 1 and value >= 0.1):
        s = '%2.1f ' % value
        snew = s#stripe10(s)
        snew = "${"+snew+"}$"
        return snew
    if(value >= 1):
        s = '%2.0f ' % value
        snew = s#stripe10(s)
        snew = "${"+snew+"}$"
        return snew
    if(value < 0.1 and value >= 0.01):
        s = '%2.2f ' % value
        snew = s#stripe10(s)
        snew = "${"+snew+"}$"
        return snew
    if(value < 0.01 and value >= 0.001):
        s = '%3.3f ' % value
        snew = s#stripe10(s)
        snew = "${"+snew+"}$"
        return snew
    else:
        s = '%3.2f ' % value
        snew = s#stripe10(s)
        snew = "${"+snew+"}$"

    return snew


def supermongolikeN5(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "${"+s+"}$"
    else:
        s = '%2d ' % value
        snew = s#stripe10(s)
        snew = "${"+snew+"}$"

    return snew

def supermongolikeN6(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "${"+s+"}$"
    else:
        s = '%2.1f ' % value
        snew = s#stripe10(s)
        snew = "${"+snew+"}$"

    return snew

def supermongolikeN7(value, pos):
    'The two args are the value and tick position'

    if(value == 0.0e0):
        s = '0'
        snew = "${"+s+"}$"
    else:
        s = '%2.2f ' % value
        snew = s#stripe10(s)
        snew = "${"+snew+"}$"

    return snew



def ScientificNotation2(s):
    # transform 1e+004 into 1e4, for example
    tup = s.split('e')
    mantissa = tup[0].rstrip('0').rstrip('.')
    sign = tup[1][0].replace('+', '')
    exponent = tup[1][1:].lstrip('0')
    if(exponent == " "):
        exponent = "0"
    if(mantissa != "1"):
        res = '%se{%s%s}' %(mantissa, sign, exponent)
        return res.replace('e{}','').replace('e',r'{\times}10^')
    else:
        res = '10^{%s%s}' %(sign,exponent)
        return res
        
def intsandonedecimal(value,pos):
    if value == int(value):
        s = '%d' % value
        snew = s
        snew = "${"+snew+"}$"
    else:
        s = '%2.1f' % value
        snew = s
        snew = "${"+snew+"}$"

    return s
