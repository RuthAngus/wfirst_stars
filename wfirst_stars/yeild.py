# Calculate the number of rotation periods we can expect to measure with
# Wfirst.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os

DATA_DIR = "data"


def teff2amp(teff):
    tb = int(teff / 500.0) * 500  # Round down to nearest 500K
    data = pd.read_csv(os.path.join(DATA_DIR, "rot_v_act{}.txt".format(tb)))
    lnps =



if __name__ == "__main__":
    teffs = np.random.randn(100)*500, + 5000
    teff2amp(4400)

    # Load the data.

# Given a teff, calculate an approximate rotation period.

# Given a teff and a period, calculate an amplitude.

# Given a magnitude, calculate the noise level.

# Based on the rotation period and amplitude, calculate whether you expect
# to be able to measure a rotation period.
