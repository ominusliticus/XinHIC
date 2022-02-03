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

##### Boiler plate plotting functions
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
###############


def PostKineticFreezeoutTransverseRadius(tau, tau_K, v_K, R_K):
    return R_K + v_K * (tau - tau_K)

def PostKineticFreezeoutVolume(tau, tau_K, v_K, R_K):
    return np.pi * PostKineticFreezeoutTransverseRadius(tau, tau_K, v_K, R_K) ** 2 * (tau + tau_K)

def XNumberDensityEvolution(tau, tau_K, v_K, R_K, init_D, init_pi, init_X, formation_rate, breakup_rate):
    
    current_X = init_X
    output = np.zeros_like(tau)
    for i, t in enumerate(tau):
        current_pi = PostKineticFreezeoutVolume(tau_K, tau_K, v_K, R_K) * init_pi / PostKineticFreezeoutVolume(t, tau_K, v_K, R_K)
        current_D = PostKineticFreezeoutVolume(tau_K, tau_K, v_K, R_K) * init_D / PostKineticFreezeoutVolume(t, tau_K, v_K, R_K)

        formation = formation_rate * (current_D ** 2) * current_pi
        breakup = breakup_rate * current_X * current_pi

        current_X += (formation - breakup) * np.diff(tau)[0]
        if current_X <= 0:
            current_X = 0

        output[i] = current_X

    return output


if __name__ == "__main__":
    gamma = 20 # MeV
    momentum = (2 * 140 * 115) ** (3. / 2.)

    formation_rate = np.array([[1, 5], [10, 50], [100, 500]]) * gamma / momentum ** 2
    breakup_rate = np.array([5, 10, 50, 100, 500]) * gamma / momentum

    initial_D = 17
    initial_pi = 700
    initial_X = 10

    R_K = 1 # fm
    v_K = np.array([0.33, 0.66, 0.99]) # units speed of light

    taus = np.linspace(0,30,300)
    
    cmap = get_cmap(10, 'tab10')
    #               densely-   loosely-dashed
    lss = ['solid', (0,(5,1)), (0,(5,10))]
    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(20,30))
    fig.patch.set_facecolor('white')
    for i in range(3): # cols of formation rate
        for j in range(2):  # rows of formation ratew
            for k in range(3):  # transvers expansion velocities
                for l in range(5):  # breakup rates
                    ax[i][j].plot(taus, XNumberDensityEvolution(taus, 1, v_K[k], R_K, initial_D, initial_pi, initial_X, formation_rate[i,j], breakup_rate[l]), color=cmap(l), ls=lss[k], lw=2)
            costumize_axis(ax[i,j], r'$\tau$ [arb units]', r'$\mathfrak n_X$ [arb units]')
            ax[i][j].set_yscale('log')

    fig.tight_layout()
    fig.savefig('proof_of_concept.pdf')
