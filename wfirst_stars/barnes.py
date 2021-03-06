import numpy as np


def tau(P, P0, t):
    x = P/P0
    return .5*(.65*t + ((.65*t)**2 - 23.4*(x**2-1) *
                        np.log(x))**.5)/np.log(x), \
        .5*(.65*t - ((.65*t)**2 - 23.4*(x**2-1)*np.log(x))**.5)/np.log(x),\



def period(age, bv):
    """
    From Barnes 2007. This is just a place holder - need to update this model.
    age in Gyr.
    Returns period in days.
    """
    return 0.7725*(bv - .4)**.601 * (age*1e3)**.5189

if __name__ == "__main__":
    print(period(10, .65))
