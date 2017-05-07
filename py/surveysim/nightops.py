"""Simulate one night of observing.
"""
from __future__ import print_function, division, absolute_import

import numpy as np

import astropy.time
import astropy.units as u

import desiutil.log

import desisurvey.exposurecalc
import desisurvey.nextobservation
import desisurvey.ephemerides
import desisurvey.config


def nightOps(night, obsplan, weather, progress, gen):
    """Simulate one night of observing.

    Use an afternoon plan, ephemerides, and simulated weather to
    schedule the observations and update the survey progress.

    Parameters
    ----------
    night : astropy.table.Row
        Row of tabulated ephemerides for this night, normally
        obtained with
        :meth:`desisurvey.ephemerides.Ephemerides.get_night`.
    obsplan : string
        Name of the file containing today's afternoon plan.
    weather : surveysim.weather.Weather
        Simulated weather conditions to use.
    progress : desisurvey.progress.Progress
        Survey progress so far, that will be updated for any
        observations taken this night.
    gen : numpy.random.RandomState
        Random number generator to use for reproducible samples.
    """
    log = desiutil.log.get_logger()
    config = desisurvey.config.Configuration()

    # Simulate the night between bright twilights.
    now = astropy.time.Time(night['brightdusk'], format='mjd')
    end_night = astropy.time.Time(night['brightdawn'], format='mjd')

    # Test if the weather permits the dome to open tonight.
    if not weather.get(now)['open']:
        log.info('Bad weather forced the dome to remain shut for the night.')
        return

    # How long to delay when we don't have a suitable target to observe.
    delay = 1. * u.day / config.num_lst_bins()

    while now < end_night:
        # Get the current weather conditions.
        conditions = weather.get(now)
        seeing, transparency = conditions['seeing'], conditions['transparency']
        # Select the next target to observe.
        target = desisurvey.nextobservation.nextFieldSelector(
            obsplan, now.mjd, progress)
        if target is None:
            # Wait until a target is available.
            now += delay
            continue
        overhead = target['overhead']
        log.debug('Selected {0} tile {1} at {2} with {3:.1f} overhead.'
                  .format(target['Program'], target['tileID'],
                          now.datetime.time(), overhead))
        # Calculate the target's airmass.
        airmass = desisurvey.utils.get_airmass(
            now, target['RA'] * u.deg, target['DEC'] * u.deg)
        # Calculate the nominal total exposure time required for this
        # target under the current observing conditions.
        total_exptime = desisurvey.exposurecalc.exposure_time(
            target['Program'], seeing, transparency, airmass, target['Ebmv'],
            night['moon_illum_frac'], target['MoonDist'],
            target['MoonAlt'])
        # Scale exposure time by the remaining SNR**2 needed for this target.
        tile = progress.get_tile(target['tileID'])
        target_exptime = total_exptime * (1 - tile['snr2frac'].sum())
        # Is this target worth observing now?
        if target_exptime > config.max_exposure_length():
            log.debug('Target {0} requires {1:.1f} exposure.  Waiting...'
                      .format(target['tileID'], target_exptime))
            now += delay
            continue
        # Calculate the number of exposures needed for cosmic ray splits.
        nexp = int(np.ceil(
            (target_exptime / config.cosmic_ray_split()).to(1).value))
        log.debug('Target {0:.1f} (total {1:.1f}) needs {2} exposures.'
                  .format(target_exptime, total_exptime, nexp))
        # Simulate the individual exposures.
        for iexp in range(nexp):
            # Advance to the start of the next exposure.
            now += overhead
            # Add random jitter with 10% RMS to the actual exposure time to
            # account for variability in the online ETC.
            exptime = target_exptime / nexp * (1 + gen.normal(scale=0.1))
            snr2frac = exptime / total_exptime
            # Record this exposure.
            progress.add_exposure(
                target['tileID'], now, exptime, snr2frac, airmass, seeing)
            now += exptime
            # Overhead for a later exposure of this target is only readout time.
            overhead = config.readout_time()
