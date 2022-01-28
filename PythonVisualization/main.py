#!/bin/python3
from __future__ import unicode_literals

from platform import uname
if 'microsoft' in uname().release:
    import os
    os.environ['MPLCONFIGDIR'] = '/tmp/'

import sys

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.ticker as tck
import matplotlib.font_manager
from matplotlib import rc, rcParams
rcParams["mathtext.fontset"]="cm"
rcParams["text.usetex"]="True"
rcParams["text.latex.preamble"]=r"\usepackage{amsmath, amsfonts}"
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
rc('text', usetex=True)

def get_cmap(n: int, name: str):
    '''
    Returns a function that maps each index in 0, 1, 2, 3, ... n-1 to a distinct
    RGB color; the keyboard argument name must be a standard mpl colormap name
    '''
    return plt.cm.get_cmap(name, n)

def costumize_axis(ax: plt.Axes, x_title: str, y_title: str):
    ax.set_xlabel(x_title, fontsize=24)
    ax.set_ylabel(y_title, fontsize=24)
    ax.tick_params(axis='both', labelsize=18, top=True, right=True)
    ax.tick_params(axis='both', which='major', direction='in', length=8)
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.tick_params(axis='both', which='minor', direction='in', length=4, top=True, right=True)
    return ax

if __name__ == '__main__':
