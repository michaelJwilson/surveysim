"""Simulation utilities that may be used by other packages.
"""
from __future__ import print_function, division, absolute_import
import numpy as np
from astropy.table import Column


def add_calibration_exposures(exposures, arcs_per_night=3, flats_per_night=3):
    """Adds calibration exposures to a set of science exposures and
    assigns exposure IDs.

    Parameters
    ----------
    exposures : :class:`astropy.table.Table`
        A table of science exposures, produced by *e.g.*
        ``desisurvey.progress.Progress.get_exposures()``.
    arcs_per_night : :class:`int`, optional
        Add this many arc exposures per night (default 3).
    flats_per_night : :class:`int`, optional
        Add this many arc exposures per night (default 3).

    Returns
    -------
    :class:`astropy.table.Table`
        A table augmented with calibration exposures.
    """
    if 'EXPID' not in exposures:
        expidc = Column(np.arange(len(exposures), dtype=np.int32),
                        name='EXPID', description='Exposure ID')
        exposures.add_column(expidc, index=0)
    for c in ('RA', 'DEC'):
        if exposures[c].unit is None:
            exposures[c].unit = 'deg'
    output = exposures[:0].copy()
    current_night = '19730703'
    expid = 0
    for i, night in enumerate(exposures['NIGHT']):
        if night != current_night:
            current_night = night
            for j in range(flats_per_night):
                output.add_row([expid, -1, -1, 0.0, 0.0, 0.0, night,
                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                'NONE', 'flat'])
                expid += 1
            for j in range(arcs_per_night):
                output.add_row([expid, -1, -1, 0.0, 0.0, 0.0, night,
                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                'NONE', 'arc'])
                expid += 1
        output.add_row(exposures[i])
        output['EXPID'][expid] = expid
        expid += 1
    return output
