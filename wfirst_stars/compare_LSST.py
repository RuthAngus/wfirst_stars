from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
import pandas as pd


def compare_pgram(fname, yr):

    # load recovered
    true_periods, true_amps = \
            np.array(pd.read_csv("simulations/{0}/all_truths.txt"
                                 .format(fname)).T)
    pfile = "results/{0}_{1}yr_results.txt".format(fname, yr)
    ids, recovered_periods = np.genfromtxt(pfile).T

    plt.clf()
    plt.plot(true_periods, recovered_periods, "k.")
    xs = np.linspace(min(true_periods), max(true_periods), 100)
    plt.plot(xs, xs, "r--")
    plt.savefig("results/pgram_compare_{0}_{1}yr".format(fname, yr))


def compare_GP(true_periods, ids, path):

    # load recovered
    recovered_periods = np.zeros_like(ids)
    errps = np.zeros_like(ids)
    # errms = np.zeros_like(ids)
    for i in range(len(ids)):
        id = str(int(ids[i])).zfill(4)
        fname = "{0}/{1}_samples.h5".format(path, id)
        if os.path.exists(fname):
            with h5py.File(fname, "r") as f:
                samples = f["samples"][...]
            nwalkers, nsteps, ndims = np.shape(samples)
            flat = np.reshape(samples, (nwalkers * nsteps, ndims))
            mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                              zip(*np.percentile(flat, [16, 50, 84], axis=0)))
            recovered_periods[i], errps[i], errps[i] = np.exp(mcmc_result[4])

    plt.clf()
    recovered_periods[recovered_periods == -9999] = 0
    plt.errorbar(true_periods, recovered_periods, yerr=errps,
                 fmt="k.", capsize=0)
    plt.ylim(0, 80)
    xs = np.linspace(min(true_periods), max(true_periods), 100)
    plt.plot(xs, xs, "r--")
    plt.savefig("GP_compare".format(path))


if __name__ == "__main__":

    # Load truths
    fname = "l45b{0}".format(sys.argv[1])
    ids, true_periods, amps = \
        np.genfromtxt("simulations/{0}/truth.txt".format(fname),
                      skip_header=1).T
    years = [1, 5, 10]
    for yr in years:
        print("year", yr)
        compare_pgram(fname, yr)
