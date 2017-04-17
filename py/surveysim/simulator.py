"""Top-level survey simulation manager.
"""
from __future__ import print_function, division, absolute_import
import numpy as np
import os.path
from shutil import copyfile
import datetime
from astropy.time import Time
from astropy.table import Table, vstack
import astropy.io.fits as pyfits
import astropy.units as u
import astropy.time
from surveysim.weather import weatherModule
from desisurvey.ephemerides import Ephemerides
from desisurvey.afternoonplan import surveyPlan
from surveysim.nightops import obsCount, nightOps
import desiutil.log
import desisurvey.utils


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
    use_jpl : bool
        Which avoidobject to use: astropy+jplephem if True, else pyephem.
    tile_file : string or None
        Name of FITS file that specifies previously observed tiles.
        The survey will start from scratch when None.
    """
    def __init__(self, start_date, stop_date, seed=None, tilesubset=None,
                 use_jpl=False, tile_file=None):
        self.log = desiutil.log.get_logger()
        self.use_jpl = use_jpl

        # Validate date range.
        self.num_days = (stop_date - start_date).days + 1
        if self.num_days <= 0:
            raise ValueError('Expected start_date < stop_date.')
        self.start_date = start_date
        self.stop_date = stop_date

        # Tabulate sun and moon ephemerides for each night of the survey.
        self.ephem = Ephemerides(start_date, stop_date, use_cache=True)

        # Build the survey plan.
        self.sp = surveyPlan(self.ephem.start.mjd, self.ephem.stop.mjd,
                             self.ephem, tilesubset=tilesubset)

        # Initialize the survey weather conditions generator.
        self.w = weatherModule(self.ephem.start.datetime, seed)

        # Resume a simulation using previously observed tiles, if requested.
        if tile_file is not None:
            tilesObserved = Table.read(tile_file)
            start_val = len(tilesObserved)+1
            self.log.info('Survey will resume after observing {0} tiles.'
                      .format(len(tilesObserved)))
        else:
            self.log.info('Survey will start from scratch.')
            tilesObserved = Table(
                names=('TILEID', 'STATUS'), dtype=('i8', 'i4'))
            tilesObserved.meta['MJDBEGIN'] = self.ephem.start.mjd
            start_val = 0
        self.tilesObserved = tilesObserved
        self.ocnt = obsCount(start_val)

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
        else:

            # Prepare a date string YYYYMMDD to use in filenames.
            date_string = '{y:04d}{m:02d}{d:02d}'.format(
                y=date.year, m=date.month, d=date.day)

            # Each day of observing starts at local noon.
            local_noon = desisurvey.utils.local_noon_on_date(date)

            # Lookup today's ephemerides.
            today = self.ephem.get(local_noon)
            sunset = today['MJDsunset']
            assert sunset > local_noon.mjd and sunset - local_noon.mjd < 1

            # Check if we are in the full-moon engineering period.
            if today['MoonFrac'] < 0.85:

                # Simulate a normal observing night.
                ntodate = len(self.tilesObserved)
                self.w.resetDome(local_noon.datetime)
                obsplan = self.sp.afternoonPlan(
                    today, date_string, self.tilesObserved)
                self.tilesObserved = nightOps(
                    today, date_string, obsplan, self.w, self.ocnt,
                    self.tilesObserved, use_jpl=self.use_jpl)
                ntiles_tonight = len(self.tilesObserved)-ntodate
                self.tiles_todo -= ntiles_tonight
                self.log.info('Observed {0} tiles tonight, {1} remaining.'
                              .format(ntiles_tonight, self.tiles_todo))
                if (self.sp.numtiles - len(self.tilesObserved)) == 0:
                    self.survey_done = True
            else:
                self.log.info(
                    'No observing around full moon ({0:.1f}% illuminated).'
                    .format(100 * today['MoonFrac']))

        self.day_index += 1
        if self.day_index == self.num_days:
            self.survey_done = True

        return not self.survey_done
