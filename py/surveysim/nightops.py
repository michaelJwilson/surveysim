"""Simulate one night of observing.
"""
from __future__ import print_function, division, absolute_import

import numpy as np

import astropy.time
import astropy.units as u

import desiutil.log

import desisurvey.etc
import desisurvey.nextobservation
import desisurvey.ephemerides
import desisurvey.config


def nightOps(night, obsplan, weather, progress, strategy, weights, gen):
    """Simulate one night of observing.

    Use an afternoon plan, ephemerides, and simulated weather to
    schedule the observations and update the survey progress.

    Parameters
    ----------
    night : astropy.table.Row
        Row of tabulated ephemerides for this night, normally
        obtained with
        :meth:`desisurvey.ephemerides.Ephemerides.get_night`.
    obsplan : string or desisurvey.plan.Planner
        Name of the file containing today's afternoon plan or a global
        planner object to use.
    weather : surveysim.weather.Weather
        Simulated weather conditions to use.
    progress : desisurvey.progress.Progress
        Survey progress so far, that will be updated for any
        observations taken this night.
    strategy : str
        Strategy to use for scheduling tiles during each night.
    weights : array or None
        Array of per-tile weights that multiply the scores used to select
        the best tile.  All weights assumed to be one when None.
    gen : numpy.random.RandomState
        Random number generator to use for reproducible samples.

    Returns
    -------
    dict
        Dictionary of total times spent during different modes (overhead, delay, live) during the night, with units.
    """
    log = desiutil.log.get_logger()
    config = desisurvey.config.Configuration()

    # Simulate the night between bright twilights.
    now = astropy.time.Time(night['brightdusk'], format='mjd')
    end_night = astropy.time.Time(night['brightdawn'], format='mjd')

    # Initialize efficiency tracking for the night.
    totals = dict(overhead=0*u.day, delay=0*u.day, live=0*u.day, dawn=0*u.day,
                  available=(end_night - now).to(u.day))

    # Test if the weather permits the dome to open tonight.
    if not weather.get(now)['open']:
        log.info('Bad weather forced the dome to remain shut for the night.')
        totals['available'] = 0 * u.day
        return totals

    # How long to delay when we don't have a suitable target to observe.
    delay_time = 1. * u.day / config.num_lst_bins()

    # Define a helper function that updates the efficiency tracking totals,
    # raises StopIteration if we reach the end of the night or else returns
    # the updated current time.
    def advance(mode, dt, totals=totals):
        if now + dt >= end_night:
            done = True
            dt = (end_night - now).to(u.s)
        else:
            done = False
        totals[mode] += dt
        log.debug('{0} {1:8s} {2:6.2f} {3}'
                  .format(now.datetime.time(), mode, dt.to(u.min), done))
        if done:
            raise StopIteration()
        return now + dt

    try:
        while True:
            # Get the current weather conditions.
            conditions = weather.get(now)
            seeing, transparency = conditions['seeing'], conditions['transparency']
            # Select the next target to observe.
            if strategy == 'baseline':
                target = desisurvey.nextobservation.nextFieldSelector(
                    obsplan, now.mjd, progress)
            else:
                target = obsplan.next_tile(
                    now, seeing, transparency, progress, strategy, weights)
            if target is None:
                log.debug('No target available at {0}. Waiting...'
                          .format(now.datetime.time()))
                now = advance('delay', delay_time)
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
            total_exptime = desisurvey.etc.exposure_time(
                target['Program'], seeing, transparency, airmass,
                target['Ebmv'], night['moon_illum_frac'], target['MoonDist'],
                target['MoonAlt'])
            # Scale exposure time by the remaining SNR needed for this target.
            tile = progress.get_tile(target['tileID'])
            target_exptime = total_exptime * max(0, 1 - tile['snr2frac'].sum())
            # Is this target worth observing now?
            if target_exptime > config.max_exposure_length():
                log.debug('Target {0} requires {1:.1f} exposure.  Waiting...'
                          .format(target['tileID'], target_exptime))
                now = advance('delay', delay_time)
                continue
            # Calculate the number of exposures needed for cosmic ray splits.
            nexp = int(np.ceil(
                (target_exptime / config.cosmic_ray_split()).to(1).value))
            log.debug('Target {0:.1f} (total {1:.1f}) needs {2} exposures.'
                      .format(target_exptime, total_exptime, nexp))
            if nexp <= 0:
                log.info('Got nexp={0} for tile {1} at {2}.'
                         .format(nexp, tile['tileID'], now.datetime.time()))
                now = advance('delay', delay_time)
                continue
            # Simulate the individual exposures.
            for iexp in range(nexp):
                # Add random jitter with 10% RMS to the target exposure time to
                # account for variability in the online ETC.
                exptime = target_exptime / nexp * (1 + gen.normal(scale=0.1))
                # Always stop by dawn.
                if now + overhead + exptime > end_night:
                    log.info('Canceling exposure {0}/{1} too close to dawn.'
                             .format(iexp + 1, nexp))
                    advance('dawn', overhead + exptime)
                # Advance to the shutter open time.
                now = advance('overhead', overhead)
                # Record this exposure.
                snr2frac = exptime / total_exptime
                progress.add_exposure(
                    target['tileID'], now, exptime, snr2frac, airmass, seeing)
                # Advance to the shutter close time.
                now = advance('live', exptime)
                # Overhead for a later exposure is only readout time.
                overhead = config.readout_time()
    except StopIteration:
        # Reached the end of the night.
        pass

    assert abs((totals['overhead'] + totals['delay'] + totals['live'] +
                totals['dawn'] - totals['available']).to(u.day).value) < 1e-6
    return totals
