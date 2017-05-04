#! /usr/bin/python
from __future__ import print_function, division, absolute_import

from datetime import datetime, timedelta
import os
from shutil import copyfile

import numpy as np

import astropy.time
import astropy.io.fits as pyfits
from astropy.table import Table, vstack
import astropy.units as u

import desiutil.log

from desisurvey.exposurecalc import expTimeEstimator
from desisurvey.nextobservation import nextFieldSelector
import desisurvey.ephemerides
import desisurvey.config


LSTres = 1.0/144.0 # Should be the same as in afternoon planner and next field selector
MaxExpLen = 3600.0 # One hour
CRsplit = 1200.0   # 20 minutes
ReadOutTime = 120.0 # Should be the same as in next field selector


def nightOps(night, date_string, obsplan, weather, progress):
    """Simulate one night of observing.

    Use an afternoon plan, ephemerides, and simulated weather to
    schedule the observations and update the survey progress.

    Parameters
    ----------
    night : astropy.table.Row
        Row of tabulated ephemerides for this night, normally
        obtained with
        :meth:`desisurvey.ephemerides.Ephemerides.get_night`.
    date_string : string
        String of the form YYYYMMDD (unused).
    obsplan : string
        Name of the file containing today's afternoon plan.
    weather : surveysim.weather.Weather
        Simulated weather conditions to use.
    progress : desisurvey.progress.Progress
        Survey progress so far, that will be updated for any
        observations taken this night.
    """
    log = desiutil.log.get_logger()
    config = desisurvey.config.Configuration()

    # Start the night during bright twilight.
    mjd = night['brightdusk']
    time = astropy.time.Time(mjd, format='mjd')

    # Test if the weather permits the dome to open tonight.
    if not weather.get(time)['open']:
        log.info('Bad weather forced the dome to remain shut for the night.')
        return

    slew = False
    ra_prev = 1.0e99
    dec_prev = 1.0e99
    while mjd < night['brightdawn']:
        # Get the current weather conditions.
        conditions = weather.get(time)
        seeing, transparency = conditions['seeing'], conditions['transparency']
        # Select the next target to observe.
        target, overhead = nextFieldSelector(
            obsplan, mjd, conditions, progress, slew, ra_prev, dec_prev)
        if target is None:
            # Wait until a target is available.
            mjd += LSTres
            slew = False
            continue
        log.debug('Selected {0} tile {1} at MJD {2:.5f} with {3:.1f}s overhead.'
                  .format(target['Program'], target['tileID'], mjd, overhead))
        # Calculate the target's airmass.
        airmass = desisurvey.utils.get_airmass(
            time, target['RA'] * u.deg, target['DEC'] * u.deg)
        # Calculate the nominal total exposure time required for this
        # target under the current observing conditions.
        total_exptime = expTimeEstimator(
            seeing, transparency, airmass, target['Program'], target['Ebmv'],
            target['DESsn2'], night['moon_illum_frac'], target['MoonDist'],
            target['MoonAlt'])
        # Calculate the target exposure time for this observation.
        # Should account for previous exposure time of partial targets here.
        target_exptime = total_exptime
        # Is this target worth observing now?
        if target_exptime > MaxExpLen:
            log.debug('Best target requires {0:.1f}s exposure.  Waiting...'
                      .format(target_exptime))
            mjd += LSTres
            slew = False # Can slew to new target while waiting.
            continue
        # Add some random jitter to the actual exposure time. This
        # should use the same generator as the weather!
        #exptime = target_exptime + np.random.normal(0.0, 20.0)
        exptime = target_exptime
        # Determine what fraction of the SNR target we have reached with
        # this exposure.
        snrfrac = exptime / total_exptime
        # Record this exposure.
        progress.add_exposure(
            target['tileID'], mjd, exptime, snrfrac, airmass, seeing)
        assert progress.get_tile(target['tileID'])['status'] == 2
        # Add extra readout time for cosmic-ray splits, if necessary.
        overhead += ReadOutTime * np.floor(exptime / CRsplit)
        # Prepare for the next exposure.
        mjd += (overhead + exptime)/86400.0
        ##tilesObserved.add_row([target['tileID'], status])
        slew = True
        ra_prev = target['RA']
        dec_prev = target['DEC']
