# Calculating the star spot variations for the Wfirst bandpass.

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si

import os

from astropy.constants import c, h, k_B


def planck(lda, T):
    exponent = h*c/(lda*k_B*T)
    spectral_radiance = (2*h*c**2/lda**5) * 1./(np.exp(exponent.value) - 1.)
    return spectral_radiance.value


def relative_fluxes(lower, upper, starteff, spotteff, fraction):
    """
    Total relative flux in the Wfirst main, IR bandpass for 10% spot
    coverage.
    """
    flux = si.quad(planck, lower, upper, args=(starteff))[0]  # star flux
    spotflux = flux*fraction \
        + (1 - fraction)*si.quad(planck, lower, upper, args=(spotteff))[0]
    return spotflux/flux


if __name__ == "__main__":
    wavelengths = np.linspace(1e-7, 5e-6, 1000)
    star = planck(wavelengths, 5000)
    starspot = planck(wavelengths, 2000)
    plt.clf()
    plt.plot(wavelengths, star)
    plt.plot(wavelengths, starspot)
    plt.plot(wavelengths, starspot + star)

    W149_l, W149_u = .927e-6, 2e-6
    Z087_l, Z087_u = .760e-6, .977e-6
    K_l, K_u = .4e-6, .9e-6
    plt.axvline(K_l, color=".5", ls="--")  # Kepler lower
    plt.axvline(K_u, color=".5", ls="--")  # Kepler upper
    plt.axvline(Z087_l, color="CornFlowerBlue")  # Z087 lower
    plt.axvline(Z087_u, color="CornFlowerBlue")  # Z087 upper
    plt.axvline(W149_l, color="HotPink")  # W149 lower
    plt.axvline(W149_u, color="HotPink")  # W149 upper

    plt.savefig("bb_spectrum")

    # Total relative flux in the Wfirst main, IR bandpass for 10% spot
    # coverage.
    W_flux = relative_fluxes(W149_l, W149_u, 5000, 2000, .9)
    print(1 - W_flux)

    # Total relative flux in the Wfirst blue bandpass
    Z_flux = relative_fluxes(Z087_l, Z087_u, 5000, 2000, .9)
    print(1 - Z_flux)

    # Total relative flux in the Kepler bandpass
    K_flux = relative_fluxes(K_l, K_u, 5000, 2000, .9)
    print(1 - K_flux)

    print((1 - W_flux)/(1 - K_flux)*100)
    print((1 - Z_flux)/(1 - K_flux)*100)
