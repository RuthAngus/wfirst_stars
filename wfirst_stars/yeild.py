# Calculate the number of rotation periods we can expect to measure with
# Wfirst.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os

DATA_DIR = "data"


def find_nearest_index(array, value):
    return (np.abs(array - value)).argmin()


def teff2amp(teff, temp_dfs, ps):
    """
    Takes a temperature and randomly generates a period based on the period-
    temperature distributions in the Buzasi files.
    Assigns an amplitude based on the period that is generated.
    Returns period and amplitude in log10 ppm. 4 is 10 ppt or 10,000 ppm or
    1%.
    Now returns amplitude in ppt
    """
    tb = int(teff / 500.) * 500  # Round down to nearest 500K
    tb = int((tb - 3500)/500)
    d = temp_dfs[tb]  # select the temperature catalogue.
    periods = ps[tb]  # select the period distribution.
    log10_period = np.random.choice(periods, 1)  # sample the period dist.
    ind = find_nearest_index(d.log10P.values, log10_period)  # select amp.
    return 10**log10_period, (10**d.log10R.values[ind])/1e3


def get_cat(teff):
    d = pd.read_csv(os.path.join(DATA_DIR, "rot_v_act{}.txt".format(teff)))
    periods, bins = [], list(d.log10P)
    bins.append(bins[-1] + .2)
    for i, bin in enumerate(d.Nbin.values):
        periods.append(np.random.uniform(bins[i], bins[i+1],
                                         d.Nbin.values[i]))
    periods = np.array([i for j in periods for i in j])
    return d, periods


def noise_calc(Hs):
    """
    Noise in ppt.
    """
    m1 =  Hs < 14.
    m2 = (14 < Hs) * (Hs < 19.6)
    m3 = Hs < 21.6
    noises = np.zeros((len(Hs)))
    noises[m1], noises[m2], noises[m3] = .6, 4, 12
    return noises


def format(df):
    """
    Reformat the dataframe data.
    """

    # select stars in the right temperature range.
    m = (3500 < df.Typ.values) * (df.Typ.values < 6500)
    truncate = 10000
    teffs = df.Typ.values[m]
    Hs = df.J.values[m]



if __name__ == "__main__":
    # Load the data.
    print("loading the main catalogue...")
    df = pd.read_csv("data.csv", index_col=None)

    # calculate the period and amplitude.
    print("loading temperature catalogues...")
    d35, p35 = get_cat(3500)
    d40, p40 = get_cat(4000)
    d45, p45 = get_cat(4500)
    d50, p50 = get_cat(5000)
    d55, p55 = get_cat(5500)
    d60, p60 = get_cat(6000)

    print("calculating teffs...")
    period, amp = [np.zeros((truncate)) for i in range(2)]
    for i, teff in enumerate(teffs):
        period[i], amp[i] = teff2amp(teff, [d35, d40, d45, d50, d55, d60],
                                     [p35, p40, p45, p50, p55, p60])

    # Given a magnitude, calculate the noise level.
    print("calculating noise...")
    noise = noise_calc(Hs)

    # Based on the rotation period, amplitude and noise level, calculate
    # whether you expect to be able to measure a rotation period.

    print("plotting...")
    plt.clf()
    plt.hist(period)
    plt.savefig("period_hist")

    plt.clf()
    plt.hist(amp)
    plt.savefig("amp_hist")

    plt.clf()
    plt.hist(noise)
    plt.savefig("noise_hist")

    m = amp > noise
    print(len(amp[m]))
