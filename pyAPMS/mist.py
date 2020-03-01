# !/usr/bin/env python3

import sys
import os
import pickle

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

def pca_extract(X, plot=True):
    """
    calculate composite score by take pca first dimension
     Args:
        Nx3 array where colums are abundance, reproducibility, and specificity
     Returns:
        N x1 array where the composite score is calculated
    """
    sc = StandardScaler()
    X_sc = sc.fit_transform(X)
    cov_mat = np.cov(X_sc.T)
    eigen_vals, eigen_vecs = np.linalg.eig(cov_mat)
    tot = sum(eigen_vals)
    var_exp = [(i/tot) for i in sorted(eigen_vals, reverse=True)]
    cum_var_exp = np.cumsum(var_exp)
    if plot:
        plt.bar(range(1, X.shape[0]+1), var_exp, alpha=0.5, align='center', label='individual explained variance')
        plt.step(range(1, X.shape[0]+1), cum_var_exp, where='mid', label='cumulative explained variance')
        plt.ylabel('Explained variance ratio')
        plt.xlabel('Principal component index')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.show()
    eigen_pairs = [(np.abs(eigen_vals[i]), eigen_vecs[:, i]) for i in range(len(eigen_vals))]
    eigen_pairs.sort(key=lambda k:k[0], reverse=True)
    # now we got the transformation matrix
    w = np.hstack((eigen_pairs[0][1][:, np.newaxis],
                eigen_pairs[1][1][:, np.newaxis]))
    return X_sc.dot(w)




def test_():
    distr = np.random.rand(100,100)
    mist = pca_extract(distr)
    print(mist)

if __name__ == '__main__':
    test_()
