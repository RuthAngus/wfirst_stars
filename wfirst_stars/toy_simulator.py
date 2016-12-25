# Generating the light curves.

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

import os

import mklc
from LSSToy import generate_visits

plotpar = {'axes.labelsize': 20,
           'xtick.labelsize': 16,
           'ytick.labelsize': 16,
           'text.usetex': True}
plt.rcParams.update(plotpar)


def simulate_LSST(id, p, a, path, noise, tmin=3, tmax=30, dur=10, plot=False):
    """ Photometry with precision of 10 ppm (?).
    Uneven time sampling that ranges (uniformily) from 3 to 30 days (?).
    Lasting 10 years (?).
    PARAMS:
    id: (str or int)
        id of the star.
    p: (float)
        Rotation period in seconds.
    a: (int)
        Amplitude in ppm
    path: (str)
        Path to save files
    tmin, tmax: (int)
        Min and max intervals between observations in days.
    dur: (int)
        Duration in years.
    noise: (int)
        Noise level (ppm). Default is 10 ppm.
    """
    print(id)
    id = str(int(id)).zfill(4)

    # The time array
    x = np.cumsum(np.random.uniform(tmin, tmax, 1000))
    x = x[x < dur * 365.25]
    x = generate_visits()
    x += -x[0]

    # a more realistic cadence model goes here.
    # x, depth = np.genfromtxt("{0}_cadence.txt".format(path)).T

    sin2incl = np.random.uniform(np.sin(0)**2, np.sin(np.pi/2)**2)
    incl = np.arcsin(sin2incl**.5)
    # tau = np.exp(np.random.uniform(np.log(p), np.log(10*p)))
    res0, res1 = mklc.mklc(x, incl=incl, p=p)
    # res0, res1 = mklc.mklc(x, p=p)

    nspot, ff, amp_err = res0
    time, area_tot, dF_tot, dF_tot0 = res1
    noise_free_y = dF_tot0 / np.median(dF_tot0) - 1
    y = noise_free_y - 1 + noise*1e-6 * np.random.randn(len(x))
    yerr = np.ones_like(y) * noise * 1e-6

    data = np.vstack((x, y, yerr))
    if not os.path.exists(path):
        os.makedirs(path)
    np.savetxt(os.path.join(path, "{}.txt".format(id)), data.T)
    truths = np.array([p, a])
    truths = np.vstack((np.array([p]), np.array([a])))
    f = open(os.path.join(path, "all_truths.txt"), "a")
    f.write("{0},{1}\n".format(truths[0][0], truths[1][0]))
    # np.savetxt(os.path.join(path, "all_truths.txt"), truths.T)
    f.close()

    if plot:
        print("plotting light curve")
        plt.clf()
        plt.errorbar(x/365.25, y, yerr=yerr, fmt="k.", capsize=0, ecolor=".5")
        plt.xlabel("$\mathrm{Time~(years)}$")
        plt.ylabel("$\mathrm{Normalised~Flux}$")
        plt.xlim(min(x/365.25), max(x/365.25))
        plt.subplots_adjust(left=.2, bottom=.12)
        plt.savefig("{0}/{1}.pdf".format(path, id))


if __name__ == "__main__":
    path = "simulations"  # where to save

    # Arrays of random (log-normal) periods and (uniform) amplitudes.
    N = 10
    ps = np.exp(np.random.uniform(np.log(2), np.log(100), N))
    amps = np.random.uniform(10, 300, N)  # ppm
    [simulate_LSST(i, ps[i], amps[i], path) for i in range(N)]

    # save the true values
    ids = np.arange(N)
    data = np.vstack((ids, ps, amps))
    np.savetxt("{0}/truth.txt".format(path), data.T)
