"""Top-level survey simulation manager.
"""
from __future__ import print_function, division, absolute_import

import datetime

import numpy as np

import astropy.table
import astropy.time
import astropy.units as u

import desiutil.log

import desisurvey.ephemerides
import desisurvey.afternoonplan
import desisurvey.plan
import desisurvey.utils
import desisurvey.config

import surveysim.nightops
import surveysim.weather


class Simulator(object):
    """Initialize a survey simulation.

    Parameters
    ----------
    start_date : datetime.date
        Survey starts on the evening of this date.
    stop_date : datetime.date
        Survey stops on the morning of this date.
    progress : desisurvey.progress.Progress
        Progress of survey at the start of this simulation.
    strategy : str
        Strategy to use for scheduling tiles during each night.
    plan : str or None
        Name of plan file to use. Required unless strategy is 'baseline'.
    seed : int or None
        Random number seed used to generate weather conditions.
    """
    def __init__(self, start_date, stop_date, progress, strategy='baseline',
                 plan=None, seed=20190823, computeHA=False):
        self.log = desiutil.log.get_logger()

        config = desisurvey.config.Configuration()
        nominal_stop_date = config.last_day()

        # Validate date range.
        self.num_days = (stop_date - start_date).days
        if self.num_days <= 0:
            raise ValueError('Expected start_date < stop_date.')
        self.start_date = start_date
        self.stop_date = stop_date

        # Tabulate sun and moon ephemerides for each night of the survey.
        if computeHA and stop_date < nominal_stop_date:
            self.ephem = desisurvey.ephemerides.Ephemerides(
                start_date, nominal_stop_date, use_cache=True)
        else:
            self.ephem = desisurvey.ephemerides.Ephemerides(
                start_date, stop_date, use_cache=True)

        if strategy == 'baseline':
            # Build the survey plan.
            self.sp = desisurvey.afternoonplan.surveyPlan(
                self.ephem.start.mjd, self.ephem.stop.mjd, self.ephem,
                computeHA)
        else:
            # Load the survey planner.
            self.sp = desisurvey.plan.Planner()
        self.strategy = strategy

        # Initialize the random number generator to use for simulating
        # the weather and adding jitter to exposure times.
        self.gen = np.random.RandomState(seed)

        # Initialize the survey weather conditions generator.
        self.weather = surveysim.weather.Weather(
            start_date, stop_date, gen=self.gen)

        if plan is not None:
            # Load the plan to use.
            config = desisurvey.config.Configuration()
            self.plan = astropy.table.Table.read(config.get_path(plan))
            assert np.all(self.sp.tiles['tileid'] == self.plan['tileid'])
        else:
            self.plan = None

        # Initialize a table for efficiency tracking.
        self.etrack = astropy.table.Table()
        self.etrack['available'] = np.zeros(self.num_days)
        self.etrack['overhead'] = np.zeros(self.num_days)
        self.etrack['delay'] = np.zeros(self.num_days)
        self.etrack['dawn'] = np.zeros(self.num_days)
        self.etrack['live'] = np.zeros(self.num_days)

        self.day_index = 0
        self.survey_done = False
        self.completed = progress.completed()
        self.progress = progress
        self.log.info(
            'Will simulate {0} to {1} with {2:.1f} / {3} tiles completed.'
            .format(start_date, stop_date, self.completed, progress.num_tiles))

    def next_day(self):
        """Simulate the next day of survey operations.

        A day runs from local noon to local noon. A survey ends, with
        this method returning False, when either we reach the last
        scheduled day or else we run out of tiles to observe.

        Returns
        -------
        bool
            True if there are more days to simulate.
        """
        if self.day_index >= self.num_days or self.survey_done:
            return False

        date = self.start_date + datetime.timedelta(days=self.day_index)
        self.log.info('Simulating {0}'.format(date))

        if desisurvey.utils.is_monsoon(date):
            self.log.info('No observing during monsoon.')
        elif self.ephem.is_full_moon(date):
            self.log.info('No observing during full moon.')
        else:

            # Each day of observing starts at local noon.
            local_noon = desisurvey.utils.local_noon_on_date(date)

            # Lookup tonight's ephemerides.
            night = self.ephem.get_night(date)

            if self.strategy == 'baseline':
                # Create the afternoon plan.
                obsplan = self.sp.afternoonPlan(night, self.progress)
            else:
                # Use the global planner.
                obsplan = self.sp

            # Simulate tonight's observing.
            totals = surveysim.nightops.nightOps(
                night, obsplan, self.weather, self.progress, self.strategy,
                self.plan, self.gen)

            # Update our efficiency tracker.
            for mode in totals:
                self.etrack[mode][self.day_index] = totals[mode].to(u.day).value

            completed = self.progress.completed()
            self.log.info(
                'Completed {0:.1f} tiles tonight, {1:.1f} remaining.'
                .format(completed - self.completed,
                        self.progress.num_tiles - completed))
            self.completed = completed
            if self.progress.num_tiles - completed < 0.1:
                self.survey_done = True

        self.day_index += 1
        if self.day_index == self.num_days:
            self.survey_done = True

        return not self.survey_done
