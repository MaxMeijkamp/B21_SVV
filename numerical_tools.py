######################################
# Authors: Nick Baquer, Paco Frantzen, Maximilian Meijkamp, Dario van Wagensveld, Mustafa Wahid and Andreas Zafiropoulos
# TU Delft, March 2020, AE3212-II Simulation, Verification and Validation group A10
# Python Script. numericalmodel.py. retrieved from https://github.com/MaxMeijkamp/A10_SVV/blob/master/numericaltools.py
# on 15 March 2020, 15:34.
######################################



import numpy as np
from functools import partial


def spline(x, f, n):
    # Spline function,
    if len(x) != len(f):
        raise ValueError("The lists are not of the same shape")
    sp_start = []
    sp_slope = []
    for i in range(0, n-1):
        splinestart = f[i]
        splineslope = ((f[i+1] - f[i])/(x[i+1] - x[i]))
        sp_start.append(splinestart)
        sp_slope.append(splineslope)
    # sp = np.vstack((np.array(sp_start), np.array(sp_slope)))
    return sp_start, sp_slope  # sp.T


def interpolate(x, f, x_target):
    # Interpolation function which gives f_target for a given x_target and x-f spline
    n = len(x)
    x = np.asarray(x)
    f = np.asarray(f)
    if x[0] > x_target or x_target > x[n-1]:
        raise ValueError("The target location is not in the range of provided function values, x_tar =", x_target,
                         "; x[0] =", x[0], "; x[n-1] =", x[n-1])
    elif x_target in x:
        return f[np.where(x_target==x)]
    else:
        sp_start, sp_slope = spline(x, f, n)
        left_i = n-2
        for i in range(n):
            if x[i] > x_target:
                left_i = i-1
                break
        f_target = sp_slope[left_i] * (x_target - x[left_i]) + sp_start[left_i]
        return f_target


def cont_spline(x_discrete, f_discrete):
    return np.vectorize(partial(interpolate, x_discrete, f_discrete))