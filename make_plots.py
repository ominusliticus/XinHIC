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
    with open('D-self-energies.dat','r') as f:
        data = np.array([[float(entry) for entry in line.split()] for line in f.readlines()])

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6 * 2,6))
    fig.patch.set_facecolor('white')
    
    cmap = ['blue', 'red'] # get_cmap(10, 'Set1')
    temps = [150, 125, 100, 75, 50]
    for i, val in enumerate([1, 3]):
        ax[0].plot(data[:,0], data[:,val], lw=3, color=cmap[i], label=r'$T=\ $' + f'{temps[i]} MeV')
        ax[1].plot(data[:,0], data[:,val+1], lw=3, color=cmap[i])
    
#        ax[0].legend(loc='lower left', fontsize=20)

#        ax[0].text(0,-750, r'$\mathfrak{n}_\pi = 2.3\times10^5$ MeV$^3$', fontsize=20)

    costumize_axis(ax[0], r'$E$ [MeV]', r'$\mathrm{Re}[\Pi_{0,\pm}(p)]$ [MeV]')
    costumize_axis(ax[1], r'$E$ [MeV]', r'$\mathrm{Im}[\Pi_{0,\pm}(p)]$ [MeV]') 
    
    fig.tight_layout()
    fig.savefig('D-self-energies.pdf')


    with open('Dt-self-energies.dat','r') as f:
        data = np.array([[float(entry) for entry in line.split()] for line in f.readlines()])

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6 * 2, 6))
    fig.patch.set_facecolor('white')

    for i, val in enumerate([1, 3]):
        ax[0].plot(data[:,0], data[:,val], lw=3, color=cmap[i], label=r'$T=\ $' + f'{temps[i]} MeV')
        ax[0].plot(data[:,0], data[:,val+1], lw=3, color=cmap[i], ls='dashed')
        ax[1].plot(data[:,0], data[:,val], lw=3, color=cmap[i], label=r'$T=\ $' + f'{temps[i]} MeV')
        ax[1].plot(data[:,0], data[:,val+1], lw=3, color=cmap[i], ls='dashed')
        
    costumize_axis(ax[0], r'$E$ [MeV]', r'$\mathrm{Re}[\Pi_{*0,\pm}(p)]$ [MeV]')
    costumize_axis(ax[1], r'$E$ [MeV]', r'$\mathrm{Im}[\Pi_{*0,\pm}(p)]$ [MeV]') 
    
    fig.tight_layout()
    fig.savefig('Dt-self-energies.pdf')
