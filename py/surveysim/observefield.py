#! /usr/bin/python

import numpy as np

# Status = 0, 1, 2 == not done, incomplete, complete
def observeField(target, exposure):

    status = 2
    real_exposure = exposure + np.random.normal(0.0, 20.0)
    realSN2 = target['DESsn2'] + np.random.uniform(0.0, 1.0)

    return status, real_exposure, realSN2

def setup_time(previous_ra, previous_dec, ra, dec):
    focus_time = 30.0
    dra = np.abs(previous_ra-ra)
    ddec = np.abs(previous_dec-dec)
    d = np.maximum(dra, ddec)
    slew_time = 11.5 + d/0.45
    overhead = focus_time + slew_time
    if overhead < 120.0:
        overhead = 120.0
    return overhead/86400.0
