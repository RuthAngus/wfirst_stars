# Generate an LSST - like light curve.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import pyfits

import os

params = {'axes.labelsize': 18,
          'font.size': 10,
          'legend.fontsize': 15,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': True}
plt.rcParams.update(params)


def load_lc(path, fnumber):
    files = sorted(glob.glob(os.path.join(path, "*slc.fits")))
    print(len(files), "files found")
    hdulist = pyfits.open(files[fnumber])
    t = hdulist[1].data
    time = t["TIME"]
    flux = t["PDCSAP_FLUX"]
    flux_err = t["PDCSAP_FLUX_ERR"]
    q = t["SAP_QUALITY"]
    m = np.isfinite(time) * np.isfinite(flux) * np.isfinite(flux_err) * \
            (q == 0)
    x = time[m]
    med = np.median(flux[m])
    y = flux[m]/med - 1
    yerr = flux_err[m]/med
    return x, y, yerr


def wfirst_cadence(x, y, yerr):
    """
    for W149, truncate to 72 days, downsample to one in 15, interpolate
    across gaps and increase white noise to 1e-3.
    Reset variability amplitude to 95%
    For Z087, same except downsample to one every 12 hours and keep
    variability amplitude..
    """

    # Truncate by keeping the end of the lc.
    x -= x[0]
    m = x < 72
    x, y, yerr = x[m], y[m], yerr[m]
    y += np.random.randn(len(y))*1e-3
    yerr = np.ones_like(yerr)*1e-3
    x1, y1, yerr1 = x[::15], y[::15], yerr[::15]
    x2, y2, yerr2 = x[::720], y[::720], yerr[::720]

    return x1, y1, yerr1, x2, y2, yerr2


def make_plot(path):
    # x0, y0, yerr0 = load_lc(path, 0)
    x1, y1, yerr1 = load_lc(path, 1)
    x2, y2, yerr2 = load_lc(path, 2)
    x3, y3, yerr3 = load_lc(path, 3)
    x = np.concatenate((x1, x2, x3))
    y = np.concatenate((y1, y2, y3))
    yerr = np.concatenate((yerr1, yerr2, yerr3))

    xW, yW, yerrW, xZ, yZ, yerrZ = wfirst_cadence(x, y, yerr)
    plt.clf()
    plt.errorbar(xW, yW, yerr=yerrW, fmt="k.", capsize=0, ecolor=".7",
                 zorder=0)
    plt.plot(xZ, yZ, "r.", zorder=1)
    plt.xlabel("$\mathrm{Time~(Days)}$")
    plt.ylabel("$\mathrm{Normalised~Flux}$")
    plt.savefig("example_lightcurve.pdf")


if __name__ == "__main__":
    # Load Kepler light curve, quarter 1
    # path = "/Users/ruthangus/.kplr/data/lightcurves/002998253"
    # path = "/Users/ruthangus/.kplr/data/lightcurves/003223000"

    path = "/Users/ruthangus/.kplr/data/lightcurves/007668623"
    make_plot(path)
    assert 0

    data = pd.read_csv("garcia.txt")
    kids = data.KID.values
    for kid in kids:
        id = str(int(kid)).zfill(9)
        print(id)
        path = "/Users/ruthangus/.kplr/data/lightcurves/{}".format(id)
        x, y, yerr = load_lc(path)
        if len(x):
            x, y, yerr = wfirst_cadence(x, y, yerr)
            plt.clf()
            plt.errorbar(x, y, yerr=yerr, fmt="k.", capsize=0, ecolor=".7")
            plt.xlabel("$\mathrm{Time~(Days)}$")
            plt.ylabel("$\mathrm{Normalised~Flux}$")
            plt.savefig("{}".format(id))
