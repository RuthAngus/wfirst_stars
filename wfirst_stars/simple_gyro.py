import numpy as np


def tau(P, P0, t):
    x = P/P0
    return .5*(.65*t + ((.65*t)**2 - 23.4*(x**2-1) *
                        np.log(x))**.5)/np.log(x), \
        .5*(.65*t - ((.65*t)**2 - 23.4*(x**2-1)*np.log(x))**.5)/np.log(x),\



def period(age, bv):
    """
    From Angus 2015. This is just a place holder - need to update this model.
    age in Gyr.
    Returns period in days.
    """
    a, b, n = .4, .31, .55
    return a*(bv - .4)**b * (age*1e3)**n

if __name__ == "__main__":
    print(period(4.5, .65))
