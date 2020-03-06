# !/usr/bin/env python3

import sys
import os
import pickle

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans


def plot_pca(X_svd):
    fig = plt.figure(1, figsize=(4, 3))
    plt.clf()
    ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

    plt.cla()
    y_pred = KMeans(n_clusters=2, random_state=0).fit_predict(X_svd)

    # Reorder the labels to have colors matching the cluster results
    ax.scatter(X_svd[:, 0], X_svd[:, 1], X_svd[:, 2], c=y_pred, edgecolor='k')

    ax.w_xaxis.set_ticklabels([])
    ax.w_yaxis.set_ticklabels([])
    ax.w_zaxis.set_ticklabels([])
    plt.show()


def calc_abu_rep(q_df):
    """
    calculate abundance and reproducibility
    abundance is defined as mean of bait-prey over all replicates bound [0,1]
    reproducibility is the normalized entropy of the vector in repl bound [0,1]
    Args:
        - q_df:  grouped dataframe for bait across multiple replicates
    Returns:
        - abu: np.array of abundance for that bait
        - rep: np.array of rep for bait
    Raises:
    """
    # works on i.e on all preys in one go
    # is grouped by bait
    q = q_df.values
    abu = np.mean(q, axis=0).flatten()
    #abu needs to have same number of rows but 1 col
    norm_entr = lambda x: np.sum(x*np.log(x))/np.log2(x.shape[0])**-1
    #entr needs to have same number of rows but 1 col
    rep = np.apply_along_axis(norm_entr, 0, q).flatten()
    metr_df = pd.DataFrame(np.array([abu, rep]), columns=list(q_df))
    metr_df = metr_df.fillna(0, inplace=False).T
    metr_df.columns = ['abu', 'repl']
    return metr_df


def create_test_data():
    """
    create input format
    """
    rnd = np.random.uniform(0, 100, size=(9,100))
    rnd = pd.DataFrame(rnd)
    rnd.index = [1,1,1,2,2,2,3,3,3]
    rnd.columns = range(0,100)
    rnd['Bait'] = rnd.index
    return rnd


def calc_mist_matrix(df, weights='PCA'):
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
        wd
    """
    df.set_index('Bait', inplace=True)
    # convert to quantities
    q = df.apply(lambda x: x/np.sum(x), axis=1)
    metr = df.groupby(['Bait']).apply(calc_abu_rep).reset_index()
    specificity = lambda x: x/np.sum(x)
    metr['sp'] = metr.groupby(['level_1'])['abu'].apply(specificity)
    X = metr[['abu','sp', 'repl']].values
    if weights == 'fixed':
        metr['cl'] = KMeans(n_clusters=2, random_state=0).fit_predict(X)
        metr['mist_score'] = metr['abu']*metr['repl']*metr['sp']
    elif weights == 'PCA':
        X = scale(X, axis=0, with_mean=True, with_std=True)
        pca = PCA()
        pca_out = pca.fit_transform(X)
        plot_pca(pca_out)
        metr['cl'] = KMeans(n_clusters=2, random_state=0).fit_predict(pca_out)
        metr['mist_score'] = pca_out[:, 0]
    return metr


def run_test():
    calc_mist_matrix(create_test_data())

if __name__ == '__main__':
    run_test()
