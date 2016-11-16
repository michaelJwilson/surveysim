#! /usr/bin/python

import numpy as np

def observeField(target, exposure):
    """
    Simulates the actual exposure and returns randum fluctuations around
    estimated exposure time and desired S/N.  Always returns a completed
    observation.

    Args:
       target: dictionnary containing 'DESsn2'
       exposure: float, estimated exposure time (seconds)

    Returns:
       status: integer, 0, 1, 2 == not done, incomplete, complete
       real_exposure: float (seconds)
       realSN2: float
    """

    status = 2
    real_exposure = exposure + np.random.normal(0.0, 20.0)
    realSN2 = target['DESsn2'] + np.random.uniform(0.0, 1.0)

    return status, real_exposure, realSN2

def setup_time(previous_ra, previous_dec, ra, dec):
    """
    Computes setup time: slew and focus (assumes readout can proceed during
    slew.

    Args:
        previous_ra: float (degrees)
        previous_dec: float (degrees)
        ra: float (degrees)
        dec: float (degrees)

    Returns:
        float, total setup time (days)
    """

    focus_time = 30.0
    dra = np.abs(previous_ra-ra)
    ddec = np.abs(previous_dec-dec)
    d = np.maximum(dra, ddec)
    slew_time = 11.5 + d/0.45
    overhead = focus_time + slew_time
    if overhead < 120.0:
        overhead = 120.0
    return overhead/86400.0

