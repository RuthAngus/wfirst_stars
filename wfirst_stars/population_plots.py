# load the stellar population data files and plot histograms of the stellar
# types

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

import os


def count_spectral_type(df):
    data = np.array(df)
    stype, weight = df.CL.values, data[:, -1]

    # F, G, K, M
    types = [5., 6., 7., 8., 9.]
    ntypes = 4
    counts = np.zeros((ntypes))
    for i in range(ntypes):
        m = (types[i] < stype) * (stype < types[i+1])
        counts[i] = sum(weight[m])
    return counts


def load_file(fname):
    df = pd.read_csv(fname, index_col=False, delim_whitespace=True,
                     skiprows=8)
    return df


def count_stars(files):
    """ how many stars of each spectral type? """
    type_count = np.zeros((len(files), 4))  # 4 spectral types
    table = load_file(files[0])
    for i, file in enumerate(files):  # for each field of view
        print(i, "of", len(files))
        df = load_file(file)
        type_count[i, :] = count_spectral_type(df)
        if i > 0:
            table = pd.concat((table, df))
    print(np.sum(type_count, axis=0))
    table.to_csv("data.csv")
    return table


def histograms():
    df = pd.read_csv("data.csv")
    plt.clf()
    plt.hist(df.Typ)
    plt.xlabel("Teff")
    plt.savefig("teff_hist.pdf")


if __name__ == "__main__":
    p1 = "/Users/ruthangus/projects/wfirst_stars/wfirst_stars/data"
    path = os.path.join(p1, "besanconWeighted/*.weighted")
    files = glob.glob(path)
    count_stars(files)
    # histograms(df)
