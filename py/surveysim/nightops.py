"""Simulate one night of observing.
"""
from __future__ import print_function, division, absolute_import

import numpy as np

import desisurvey.utils
import desisurvey.etc
import desisurvey.plots


def simulate_night(night, scheduler, stats, explist, weather,
                   use_twilight=False, plot=False, verbose=False):
    """Replaces nightOpsDeprecated().
    """
    night = desisurvey.utils.get_date(night)
    nightstats = stats.get_night(night)
    label = str(night)
    if verbose: print('Simulating {}.'.format(night))
    # Lookup this night's sunset and sunrise MJD values.
    night_ephem = scheduler.ephem.get_night(night)
    if use_twilight:
        begin = night_ephem['brightdusk']
        end = night_ephem['brightdawn']
    else:
        begin = night_ephem['dusk']
        end = night_ephem['dawn']
    nightstats['tsched'] = end - begin
    if verbose:
        print('begin: {}, end: {}'.format(begin, end))

    # Find weather time steps that cover this night.
    weather_mjd = weather._table['mjd'].data
    ilo = np.searchsorted(weather_mjd, begin, side='left') - 1
    ihi = np.searchsorted(weather_mjd, end, side='right') + 2
    assert weather_mjd[ilo] < begin and weather_mjd[ihi] > end
    weather_mjd = weather_mjd[ilo:ihi]
    seeing = weather._table['seeing'].data[ilo:ihi]
    transp = weather._table['transparency'].data[ilo:ihi]
    # Fix this in the weather generator instead?
    transp = np.maximum(0.1, transp)
    dome = weather._table['open'].data[ilo:ihi]

    if not np.any(dome):
        if verbose: print('Dome closed all night.')
        return

    scheduler.init_night(night, use_twilight=use_twilight, verbose=verbose)
    ETC = desisurvey.etc.ExposureTimeCalculator(save_history=plot)
    nexp_last = explist.nexp

    # Build linear interpolators for observing conditions.
    # This implementation is faster than scipy.interpolate.interp1d()
    # when mjd values are gradually increasing.
    weather_idx = 0
    dmjd_weather = weather_mjd[1] - weather_mjd[0]
    def get_weather(mjd):
        nonlocal weather_idx
        while mjd >= weather_mjd[weather_idx + 1]:
            weather_idx += 1
        s = (mjd - weather_mjd[weather_idx]) / dmjd_weather
        return (
            seeing[weather_idx] * (1 - s) + seeing[weather_idx + 1] * s,
            transp[weather_idx] * (1 - s) + transp[weather_idx + 1] * s)

    # Define time intervals to use in units of days (move to config?)
    NO_TILE_AVAIL_DELAY = 30. / 86400.

    # Step through the night.
    dome_is_open = False
    mjd_now = weather_mjd[0]
    completed_last = scheduler.completed_by_pass.copy()
    while mjd_now < end:
        if not dome_is_open:
            # Advance to the next dome opening, if any.
            idx_now = np.searchsorted(weather_mjd, mjd_now, side='left')
            if not np.any(dome[idx_now:]):
                # Dome is closed for the rest of the night.
                mjd_now = end
                break
            idx_open = idx_now + np.argmax(dome[idx_now:])
            assert dome[idx_open] == True and (idx_open == 0 or dome[idx_open - 1] == False)
            mjd_now = weather_mjd[idx_open]
            if mjd_now >= end:
                # The next dome opening is after the end of the night.
                # This can happen if we are not using twilight.
                break
            # Find the next closing.
            if np.all(dome[idx_open:]):
                next_dome_closing = end
            else:
                idx_close = idx_open + np.argmin(dome[idx_open:])
                assert dome[idx_close] == False and dome[idx_close - 1] == True
                next_dome_closing = min(end, weather_mjd[idx_close])
            dome_is_open = True
            weather_idx = idx_open

        # == NEXT TILE ===========================================================        
        # Dome is open from mjd_now to next_dome_closing.
        mjd_last = mjd_now
        tdead = 0.
        # Get the current observing conditions.
        seeing_now, transp_now = get_weather(mjd_now)
        # Get the next tile to observe from the scheduler.
        tileid, passnum, snr2frac_start, exposure_factor, program, mjd_program_end = scheduler.next_tile(
            mjd_now, ETC, seeing_now, transp_now)
        if tileid is None:
            # Deadtime while we delay and try again.
            mjd_now += NO_TILE_AVAIL_DELAY
            if mjd_now >= next_dome_closing:
                # Dome closed during deadtime.
                mjd_now = next_dome_closing
                dome_is_open = False
            tdead += mjd_now - mjd_last
        else:
            # Setup for a new field.
            mjd_now += ETC.NEW_FIELD_SETUP
            if mjd_now >= next_dome_closing:
                # Setup interrupted by dome closing.
                mjd_now = next_dome_closing
                dome_is_open = False
                # Record an aborted setup.
                nightstats['nsetup_abort'][passnum] += 1
            else:
                # Record a completed setup.
                nightstats['nsetup'][passnum] += 1
            # Charge this as setup time whether or not it was aborted.
            nightstats['tsetup'][passnum] += mjd_now - mjd_last

            if dome_is_open:            
                # Lookup the program name for this pass.
                program = scheduler.pass_program[passnum]
                # Loop over repeated exposures of the same tile.
                continue_this_tile = True
                while continue_this_tile:
                    # -- NEXT EXPOSURE ---------------------------------------------------
                    # Get the current observing conditions.
                    seeing_now, transp_now = get_weather(mjd_now)
                    sky_now = 1.
                    # Use the ETC to control the shutter.
                    mjd_open_shutter = mjd_now
                    ETC.start(mjd_now, tileid, program, snr2frac_start, exposure_factor,
                              seeing_now, transp_now, sky_now)
                    integrating = True
                    while integrating:
                        mjd_now += ETC.UPDATE_INTERVAL
                        if mjd_now >= next_dome_closing:
                            # Current exposure is interrupted by dome closing.
                            mjd_now = next_dome_closing
                            dome_is_open = False
                            integrating = False
                            continue_this_tile = False
                        elif mjd_now >= mjd_program_end:
                            # Current exposure is interrupted by a program change.
                            mjd_now = mjd_program_end
                            integrating = False
                            continue_this_tile = False
                        # Get the current observing conditions.
                        seeing_now, transp_now = get_weather(mjd_now)
                        sky_now = 1.
                        # Update the SNR.
                        if not ETC.update(mjd_now, seeing_now, transp_now, sky_now):
                            # Current exposure reached its target SNR according to the ETC.
                            integrating= False
                    # stop() will return False if this is a cosmic split and
                    # more integration is still required.
                    if ETC.stop(mjd_now):
                        continue_this_tile = False

                    # Record this exposure
                    assert np.allclose(ETC.exptime, mjd_now - mjd_open_shutter)
                    nightstats['tscience'][passnum] += ETC.exptime
                    nightstats['nexp'][passnum] += 1
                    explist.add(
                        mjd_now - ETC.exptime, 86400 * ETC.exptime,
                        tileid, passnum,
                        snr2frac_start, ETC.snr2frac, seeing_now, transp_now, sky_now)
                    scheduler.update_tile(tileid, ETC.snr2frac)

                    if continue_this_tile:
                        # Prepare for the next exposure of the same tile.
                        snr2frac_start = ETC.snr2frac
                        mjd_split_start = mjd_now
                        mjd_now += ETC.SAME_FIELD_SETUP
                        if mjd_now >= next_dome_closing:
                            # Setup for next exposure of same tile interrupted by dome closing.
                            mjd_now = next_dome_closing
                            dome_is_open = False
                            continue_this_tile = False
                            # Record an aborted split.
                            nightstats['nsplit_abort'][passnum] += 1
                        else:
                            # Record a completed split.
                            nightstats['nsplit'][passnum] += 1
                        # Charge this as split time, whether or not is was aborted.
                        nightstats['tsplit'][passnum] += mjd_now - mjd_split_start
                    # --------------------------------------------------------------------

        # Update statistics for this program.
        pidx = scheduler.program_index[program]
        nightstats['tdead'][pidx] += tdead
        nightstats['topen'][pidx] += mjd_now - mjd_last

        # All done if we have observed all tiles.
        if scheduler.complete():
            break
        # ========================================================================

    # Save the number of tiles completed per pass in the nightly statistics.
    nightstats['completed'][:] = scheduler.completed_by_pass - completed_last

    if plot:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
        ax = axes[0]
        ax.plot(weather_mjd, seeing, 'r-', label='Seeing')
        ax.plot([], [], 'b-', label='Transparency')
        ax.legend(ncol=2, loc='lower center')
        ax.set_ylabel('Seeing FWHM [arcsec]')
        rhs = ax.twinx()
        rhs.plot(weather_mjd, transp, 'b-')
        rhs.set_ylabel('Transparency')
        ax = axes[1]
        changes = np.where(np.abs(np.diff(dome)) == 1)[0]
        for idx in changes:
            ax.axvline(weather_mjd[idx + 1], ls='-', c='r')
        mjd_history = np.array(ETC.history['mjd'])
        snr2frac_history = np.array(ETC.history['snr2frac'])
        for expinfo in explist._exposures[nexp_last: explist.nexp]:
            color = desisurvey.plots.program_color[scheduler.pass_program[expinfo['passnum']]]
            t1 = expinfo['mjd']
            t2 = t1 + expinfo['exptime'] / 86400
            y1, y2 = expinfo['snr2frac_start'], expinfo['snr2frac_stop']
            sel = (mjd_history >= t1) & (mjd_history <= t2)
            ax.fill_between(mjd_history[sel], snr2frac_history[sel], color=color, alpha=0.5, lw=0)
        for t in scheduler.night_changes:
            ax.axvline(t, c='b', ls=':')
        ax.set_xlim(weather_mjd[0], weather_mjd[-1])
        ax.set_xlabel('MJD During {}'.format(label))
        ax.set_ylim(0, 1)
        ax.set_ylabel('Integrated SNR2 Fraction')


import astropy.time
import astropy.units as u

import desiutil.log

import desisurvey.config


def nightOpsDeprecated(date, ephem, scheduler, weather, progress, strategy, plan, scores,
             gen):
    """Simulate one night of observing. Superceded by simulate_night().

    Use an afternoon plan, ephemerides, and simulated weather to
    schedule the observations and update the survey progress.

    Parameters
    ----------
    date : datetime.date
        Date when this night starts.
    ephem : desisurvey.ephemerides.Ephemerides
        Tabulated ephemerides data to use for simulating this night.
    scheduler : desisurvey.old.schedule.Scheduler
        Scheduler object to use for selecting next tiles.
    weather : surveysim.weather.Weather
        Simulated weather conditions to use.
    progress : desisurvey.progress.Progress
        Survey progress so far, that will be updated for any
        observations taken this night.
    strategy : str
        Strategy to use for scheduling tiles during each night.
    plan : astropy.table.Table
        Table that specifies tile priorities and design hour angles.
    scores : list or None
        Append an array of per-tile scheduler scores to this list for each
        exposure unless None. Scores are saved as float32 values.
    gen : numpy.random.RandomState
        Random number generator to use for reproducible samples.

    Returns
    -------
    dict
        Dictionary of total times spent during different modes (overhead, delay, live) during the night, with units.
    """
    log = desiutil.log.get_logger()
    config = desisurvey.config.Configuration()

    # Lookup tonight's ephemerides.
    night = ephem.get_night(date)

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
    delay_time = 1. * u.min

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
            seeing, transparency = (
                conditions['seeing'], conditions['transparency'])
            if transparency < 0.05:
                log.warn('Clipping transparency {0:.6f} to 0.05'
                         .format(transparency))
                transparency = 0.05
            # Select the next target to observe.
            target = scheduler.next_tile(
                now, ephem, seeing, transparency, progress, strategy, plan)
            if target is None:
                log.debug('No target available at {0}. Waiting...'
                          .format(now.datetime.time()))
                now = advance('delay', delay_time)
                continue
            overhead = target['overhead']
            log.info('Selected {0} tile {1} at {2} with {3:.1f} overhead.'
                     .format(target['Program'], target['tileID'],
                              now.datetime.time(), overhead))
            # Calculate the target's airmass.
            airmass = desisurvey.utils.get_airmass(
                now, target['RA'] * u.deg, target['DEC'] * u.deg)
            # Calculate the nominal total exposure time required for this
            # target under the current observing conditions.
            moonfrac = night['moon_illum_frac']
            moonsep = target['MoonDist']
            moonalt = target['MoonAlt']
            total_exptime = desisurvey.etc.exposure_time(
                target['Program'], seeing, transparency, airmass,
                target['Ebmv'], moonfrac, moonsep, moonalt)
            # Scale exposure time by the remaining SNR needed for this target.
            tile = progress.get_tile(target['tileID'])
            target_exptime = total_exptime * max(0, 1 - tile['snr2frac'].sum())
            # Do not re-observe a target that has already been completed.
            if target_exptime == 0:
                log.info('Tile {0} already completed at {1}.'
                         .format(tile['tileid'], now.datetime.time()))
                now = advance('delay', delay_time)
                continue
            log.info('Target exposure time {0:.1f} = {1:.1f}.'
                     .format(target_exptime.to(u.s), target_exptime.to(u.min)))
            # Clip the exposure time if necessary.
            if target_exptime > config.max_exposure_length():
                log.info('Clip exposure time {0:.1f} -> {1:.1f} for tile {2}.'
                         .format(target_exptime, config.max_exposure_length(),
                                 target['tileID']))
                target_exptime = config.max_exposure_length()
            # Calculate the number of exposures needed for cosmic ray splits.
            nexp = int(np.ceil(
                (target_exptime / config.cosmic_ray_split()).to(1).value))
            log.debug('Target {0:.1f} (total {1:.1f}) needs {2} exposures.'
                      .format(target_exptime, total_exptime, nexp))
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
                    target['tileID'], now, exptime, snr2frac, airmass, seeing,
                    transparency, moonfrac, moonalt, moonsep)
                if scores is not None:
                    scores.append(target['score'].astype(np.float32))
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
