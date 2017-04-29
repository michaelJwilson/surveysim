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
import desisurvey.utils

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
    seed : int or None
        Random number seed used to generate weather conditions.
    tilesubset : array or None
        Array of integer tileIDs to use while ignoring others
        in the DESI footprint.
    tile_file : string or None
        Name of FITS file that specifies previously observed tiles.
        The survey will start from scratch when None.
    """
    def __init__(self, start_date, stop_date, seed=None, tilesubset=None,
                 tile_file=None):
        self.log = desiutil.log.get_logger()

        # Validate date range.
        self.num_days = (stop_date - start_date).days
        if self.num_days <= 0:
            raise ValueError('Expected start_date < stop_date.')
        self.start_date = start_date
        self.stop_date = stop_date

        # Tabulate sun and moon ephemerides for each night of the survey.
        self.ephem = desisurvey.ephemerides.Ephemerides(
            start_date, stop_date, use_cache=True)

        # Build the survey plan.
        self.sp = desisurvey.afternoonplan.surveyPlan(
            self.ephem.start.mjd, self.ephem.stop.mjd, self.ephem,
            tilesubset=tilesubset)

        # Initialize the survey weather conditions generator.
        self.weather = surveysim.weather.Weather(
            start_date, stop_date, seed=seed)

        # Resume a simulation using previously observed tiles, if requested.
        if tile_file is not None:
            tilesObserved = Table.read(tile_file)
            start_val = len(tilesObserved)+1
            self.log.info('Survey will resume after observing {0} tiles.'
                      .format(len(tilesObserved)))
        else:
            self.log.info('Survey will start from scratch.')
            tilesObserved = astropy.table.Table(
                names=('TILEID', 'STATUS'), dtype=('i8', 'i4'))
            tilesObserved.meta['MJDBEGIN'] = self.ephem.start.mjd
            start_val = 0
        self.tilesObserved = tilesObserved
        self.ocnt = surveysim.nightops.obsCount(start_val)

        self.day_index = 0
        self.survey_done = False
        self.tiles_todo = self.sp.numtiles
        self.log.info(
            'Simulator initialized for {0} to {1} with {2} tiles remaining.'
            .format(start_date, stop_date, self.tiles_todo))


    def next_day(self):
        """Simulate the next day of survey operations.

        A day runs from local noon to local noon. A survey ends, with this
        method returning False, when either we reach the last scheduled day or
        else we run out of tiles to observe.

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

            # Prepare a date string YYYYMMDD to use in filenames.
            date_string = '{y:04d}{m:02d}{d:02d}'.format(
                y=date.year, m=date.month, d=date.day)

            # Each day of observing starts at local noon.
            local_noon = desisurvey.utils.local_noon_on_date(date)

            # Lookup tonight's ephemerides.
            night = self.ephem.get_night(date)

            # Simulate a normal observing night.
            ntodate = len(self.tilesObserved)
            obsplan = self.sp.afternoonPlan(
                night, date_string, self.tilesObserved)
            self.tilesObserved = surveysim.nightops.nightOps(
                night, date_string, obsplan, self.weather, self.ocnt,
                self.tilesObserved)
            ntiles_tonight = len(self.tilesObserved)-ntodate
            self.tiles_todo -= ntiles_tonight
            self.log.info('Observed {0} tiles tonight, {1} remaining.'
                          .format(ntiles_tonight, self.tiles_todo))
            if (self.sp.numtiles - len(self.tilesObserved)) == 0:
                self.survey_done = True

        self.day_index += 1
        if self.day_index == self.num_days:
            self.survey_done = True

        return not self.survey_done
