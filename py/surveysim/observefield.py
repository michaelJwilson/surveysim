from __future__ import print_function, division
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
