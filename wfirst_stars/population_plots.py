# load the stellar population data files and plot histograms of the stellar
# types

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

import os

# plt.rcParams.update({'axes.labelsize': 18, 'font.size': 10,
#                      'legend.fontsize': 15, 'xtick.labelsize': 18,
#                      'ytick.labelsize': 18, 'text.usetex': True})


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


def join_stars(files):
    table = load_file(files[0])
    for i, file in enumerate(files):  # for each field of view
        print(i, "of", len(files))
        df = load_file(file)
        if i > 0:
            table = pd.concat((table, df))
    table.to_csv("data.csv")
    return table


def histograms(df):

    teffs = df.Typ.values
    distances = df.Dist.values
    logg = df.logg.values

    m = (2000 < teffs) * (teffs < 7000)
    plt.clf()
    print(len(teffs[m]))
    plt.hist(teffs[m], 10)
    plt.xlim(7000, 2000)
    # plt.xlabel("$\mathrm{T}_{\mathrm{eff}}$")
    plt.xlabel("Teff")
    plt.savefig("teff_hist.pdf")

    m = (3000 < teffs) * (teffs < 7000) * (logg < 4.)
    plt.clf()
    print(len(distances[m]))
    plt.hist(distances[m], 20, histtype="stepfilled", color="w")
    # plt.xlabel("$\mathrm{Distance~(Kpc)}$")
    plt.xlabel("Distance (Kpc)")
    plt.savefig("dist_hist.pdf")


if __name__ == "__main__":
    # p1 = "/Users/ruthangus/projects/wfirst_stars/wfirst_stars/data"
    # path = os.path.join(p1, "besanconWeighted/*.weighted")
    # df = join_stars(glob.glob(path))

    df = pd.read_csv("data.csv")
    histograms(df)
    # count_spectral_type(df)
