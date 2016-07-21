#! /usr/bin/python

import numpy as np

# Status = 0, 1, 2 == not done, incomplete, complete
def observeField(target, exposure):

    status = 2
    real_exposure = exposure + np.random.normal(0.0, 20.0)
    realSN2 = target['DESsn2'] + np.random.uniform(0.0, 1.0)

    return status, real_exposure, realSN2
