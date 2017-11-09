"""Simulation utilities that may be used by other packages.
"""
from __future__ import print_function, division, absolute_import
import numpy as np
from astropy.table import Column


def add_calibration_exposures(exposures, flats_per_night=3, arcs_per_night=3,
                              darks_per_night=0, zeroes_per_night=0,
                              exptime=None, readout=30.0):
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
    exptime : :class:`dict`, optional
        A dictionary setting calibration exposure times for each
        calibration flavor.
    readout : :class:`float`, optional
        Set readout time for calibration exposures (default 30.0 s).

    Returns
    -------
    :class:`astropy.table.Table`
        A table augmented with calibration exposures.
    """
    if exptime is None:
        exptime = {'flat': 10.0, 'arc': 10.0, 'dark': 1000.0, 'zero': 0.0}
    output = exposures[:0].copy()
    expid_in_exposures = 'EXPID' in exposures.colnames
    for c in ('RA', 'DEC'):
        if output[c].unit is None:
            output[c].unit = 'deg'
    current_night = '19730703'
    expid = 0
    calib_time = lambda x: exptime[x] + readout
    calib_sequence = (['flat']*flats_per_night + ['arc']*arcs_per_night +
                      ['dark']*darks_per_night + ['zero']*zeroes_per_night)
    calib_times = np.cumsum(np.array([calib_time(c)
                                      for c in calib_sequence]))[::-1]
    # None fields get filled in for each calibration exposure.
    calib_data = {'TILEID': -1,
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
    if expid_in_exposures:
        calib_data['EXPID'] = None
    for i, night in enumerate(exposures['NIGHT']):
        if night != current_night:
            current_night = night
            for j, c in enumerate(calib_sequence):
                if expid_in_exposures:
                    calib_data['EXPID'] = expid
                calib_data['NIGHT'] = night
                calib_data['MJD'] = exposures['MJD'][i] - calib_times[j]/86400.0
                calib_data['EXPTIME'] = exptime[c]
                calib_data['FLAVOR'] = c
                output.add_row(calib_data)
                expid += 1
        output.add_row(exposures[i])
        if expid_in_exposures:
            output['EXPID'][expid] = expid
        expid += 1
    if not expid_in_exposures:
        output.add_column(Column(np.arange(len(output), dtype=np.int32),
                                 name='EXPID', description='Exposure ID'),
                          index=0)
    return output
