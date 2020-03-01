# !/usr/bin/env python3

import sys
import os
import pickle

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def normalize_wd(wd_arr, norm_=0.98):
    """
    normalize
    Args:
        arr numpy array of wd scores
        norm_ [0,1] number corresponding to the wd_arr quantile to norm default 0.98
    Returns:
        normalized value
    """
    return wd_arr / np.quantile(wd_arr[wd_arr > 0], norm_)


def vec_wd_score(arr, norm):
    """
    vectorized wd score
    Args:
        arr: 1d array array
        norm: boolean for normalization or not
    Returns:
        a single wd score for this row
    Raises:
    """
    pres = arr[arr>0].shape[0]
    npres = arr[arr==0].shape[0]
    ntot =pres+npres
    mu_ = np.sum(arr)/ntot
    sum_sq_err = np.sum((arr - mu_)**2) + ((mu_**2) * npres)
    sd_prey = np.sqrt(sum_sq_err / (ntot-1))
    wj = sd_prey / mu_
    if wj < 1:
        wj = 1
    wd_inner = (ntot / pres) * wj
    wd = arr * wd_inner
    if norm:
        return normalize_wd(wd)
    else:
        return wd


def plot_distr(wd, sim, cutoff, plotname):
    """
    plot of real and simulated data
     Args:
        wd is matrix of wd scores from calc_wd
        sim is distribution of simulated scores
        cutoff is npquantile of sim distr
        plotname is used to assess if save to file or visualize
     Returns:
        None
    """
    wd = wd.flatten()
    sim = sim.flatten()
    plt.figure(figsize=(6, 6))
    binNum = 100
    dist = np.unique(np.concatenate((wd, sim)))
    binwidth = (max(dist) - min(dist)) / binNum
    plt.hist(wd, bins=np.arange(min(dist), max(dist) + binwidth, binwidth),
             color='r', edgecolor='r', alpha=0.3, label='Target')
    plt.hist(sim, bins=np.arange(min(dist), max(dist) + binwidth, binwidth),
             color='b', edgecolor='b', alpha=0.3, label='Simulated')
    plt.ylabel('Frequency')
    plt.axvline(x=cutoff,color='green')
    plt.legend()
    plt.tight_layout()
    if plotname:
        plt.savefig(plotname, dpi=800, bbox_inches='tight')
    else:
        plt.show()
    plt.close()
    return True


def calc_wd_matrix(m, iteration=1000, q=0.9, norm=False, plot=False):
    """
    get a NxM matrix and calculate wd then for iteration creates dummy matrix
    and score them to get distribution of simulated interactors for each bait
    Args:
        m stats matrix following http://besra.hms.harvard.edu/ipmsmsdbs/cgi-bin/tutorial.cgi format
        iteration number of iteration for generating simulated distribution
        quantile to filter interactors for
        norm quantile based normalization of interaction
        plot boolean for plotting distribution of real and simulated data
    Returns:
        wd scores matrix
    """
    wd = np.array([vec_wd_score(m[i], norm) for i in range(m.shape[1])])
    i = 0
    rand_dist = []
    p = m.flatten()
    while i <= iteration:
        np.random.shuffle(p)
        p_arr = np.random.choice(p, m.shape[0])
        # force to have some numbers inside
        while not np.any(p_arr):
            np.random.choice(p, m.shape[0])
        # print('iteration {} of {}'.format(i, iteration))
        rand_dist.append(vec_wd_score(p_arr, norm).flatten())
        i+=1
    rand_dist = np.array(rand_dist).reshape(-1,1)
    cutoff = np.quantile(rand_dist[rand_dist>0], q)
    if plot:
        plot_distr(wd, rand_dist, cutoff, None)
    wd[wd < cutoff] = 0
    return wd

def test_():
    distr = np.random.rand(100,100)
    wd_scores = calc_wd_matrix(distr, iteration=1000, q=0.8, plot=True)

if __name__ == '__main__':
    test_()
