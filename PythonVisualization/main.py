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

from parameters import Parameters
import SelfEnergies as SE

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
    params = Parameters()

    D_momentum = np.linspace(0.1, 140, 1000)
    # pion_density = np.array([0.03 * (200 ** 3), 3e-4 * (200 ** 3)])
    pion_density = np.array([0.03 * (200 ** 3)])
    Ts = np.linspace(50, 150, 10, endpoint=True)

    Dzero_ReSelfEnergy = np.zeros((Ts.size,pion_density.size,D_momentum.size))
    Dzero_ImSelfEnergy = np.zeros((Ts.size,pion_density.size,D_momentum.size))
    Dplus_ReSelfEnergy = np.zeros((Ts.size,pion_density.size,D_momentum.size))
    Dplus_ImSelfEnergy = np.zeros((Ts.size,pion_density.size,D_momentum.size))

    for n in range(Ts.size):
        params.T = Ts[n]
        for i in range(pion_density.size):
            for j in range(D_momentum.size):
                Dzero_ReSelfEnergy[n,i,j] = SE.ReDmesonSelfEnergy(D_momentum[j], params, pion_density[i], '0')
                Dzero_ImSelfEnergy[n,i,j] = SE.ImDmesonSelfEnergy(D_momentum[j], params, pion_density[i], '0')
                Dplus_ReSelfEnergy[n,i,j] = SE.ReDmesonSelfEnergy(D_momentum[j], params, pion_density[i], '+')
                Dplus_ImSelfEnergy[n,i,j] = SE.ImDmesonSelfEnergy(D_momentum[j], params, pion_density[i], '+')

    linesytles = ['solid','dashed']
    cmap = get_cmap(10, 'tab10')
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(6 * 2, 6 * 2))
    fig.patch.set_facecolor('white')
    for n in range(Ts.size):
        for i in range(pion_density.size):
            ax[0,0].plot(D_momentum, Dzero_ReSelfEnergy[n,i], lw=2, color=cmap(n), ls=linesytles[i])
            ax[1,0].plot(D_momentum, Dplus_ReSelfEnergy[n,i], lw=2, color=cmap(n), ls=linesytles[i])
            ax[0,1].plot(D_momentum, Dzero_ImSelfEnergy[n,i], lw=2, color=cmap(n), ls=linesytles[i])
            ax[1,1].plot(D_momentum, Dplus_ImSelfEnergy[n,i], lw=2, color=cmap(n), ls=linesytles[i])
    
    costumize_axis(ax[0,0], r'$p$ [MeV]', r'Re$[\Pi_0(p)]$ [MeV]')
    costumize_axis(ax[0,1], r'$p$ [MeV]', r'Im$[\Pi_0(p)]$ [MeV]')
    costumize_axis(ax[1,0], r'$p$ [MeV]', r'Re$[\Pi_+(p)]$ [MeV]')
    costumize_axis(ax[1,1], r'$p$ [MeV]', r'Im$[\Pi_+(p)]$ [MeV]')

    fig.tight_layout()
    fig.savefig('plots/D-meson-self-energy.pdf')
