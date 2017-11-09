"""Simulation utilities that may be used by other packages.
"""
from __future__ import print_function, division, absolute_import
import numpy as np
from astropy.table import Column


def add_calibration_exposures(exposures, flats_per_night=3, arcs_per_night=3,
                              darks_per_night=0, zeroes_per_night=0,
                              exptime=10.0, readout=30.0):
    """Adds calibration exposures to a set of science exposures and
    assigns exposure IDs.

    Parameters
    ----------
    exposures : :class:`astropy.table.Table`
        A table of science exposures, produced by *e.g.*
        ``desisurvey.progress.Progress.get_exposures()``.
    flats_per_night : :class:`int`, optional
        Add this many arc exposures per night (default 3).
    arcs_per_night : :class:`int`, optional
        Add this many arc exposures per night (default 3).
    darks_per_night : :class:`int`, optional
        Add this many dark exposures per night (default 0).
    zeroes_per_night : :class:`int`, optional
        Add this many zero exposures per night (default 0).
    exptime : :class:`float`, optional
        Set calibration exposure times (default 10.0 s).
    readout : :class:`float`, optional
        Set readout time for calibration exposures (default 30.0 s).

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
    calib_exptime = lambda x: 0.0 if x == 'zero' else exptime
    calib_time = lambda x: calib_exptime(x) + readout
    calib_sequence = (['flat']*flats_per_night + ['arc']*arcs_per_night +
                      ['dark']*darks_per_night + ['zero']*zeroes_per_night)
    calib_times = np.cumsum(np.array([calib_time(c)
                                      for c in calib_sequence])[::-1])
    # None fields get filled in for each calibration exposure.
    calib_data = {'EXPID': None,
                  'TILEID': -1,
                  'PASS': -1,
                  'RA': 0.0,
                  'DEC': 0.0,
                  'EBMV': 0.0,
                  'NIGHT': None,
                  'MJD': None,
                  'EXPTIME': None,
                  'SEEING': 0.0,
                  'TRANSPARENCY': 0.0,
                  'AIRMASS': 0.0,
                  'MOONFRAC': 0.0,
                  'MOONALT': 0.0,
                  'MOONSEP': 0.0,
                  'PROGRAM': 'CALIB',
                  'FLAVOR': None
                  }
    for i, night in enumerate(exposures['NIGHT']):
        if night != current_night:
            current_night = night
            for j, c in enumerate(calib_sequence):
                calib_data['EXPID'] = expid
                calib_data['NIGHT'] = night
                calib_data['MJD'] = exposures['MJD'][i] - calib_times[j]/86400.0
                calib_data['EXPTIME'] = calib_exptime(c)
                calib_data['FLAVOR'] = c
                output.add_row(calib_data)
                expid += 1
        output.add_row(exposures[i])
        output['EXPID'][expid] = expid
        expid += 1
    return output
