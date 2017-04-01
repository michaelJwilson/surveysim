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
from desisurvey.nightcal import getCalAll
from desisurvey.afternoonplan import surveyPlan
from surveysim.nightops import obsCount, nightOps
import desiutil.log


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

        # Note 1900 UTC is midday at KPNO, which is in Mountain Standard Time
        # UTC-7 and does not observe daylight savings.
        local_noon = datetime.time(hour=19)
        self.startdate = Time(datetime.datetime.combine(start_date, local_noon))
        self.enddate = Time(datetime.datetime.combine(stop_date, local_noon))
        self.log.info('Simulator initialized for {0} to {1}'.format(
            self.startdate, self.enddate))

        # Tabulate sun and moon ephemerides for each night of the survey.
        self.surveycal = getCalAll(self.startdate, self.enddate, use_cache=True)

        # Build the survey plan.
        self.sp = surveyPlan(self.startdate.mjd, self.enddate.mjd,
                             self.surveycal, tilesubset=tilesubset)

        # Initialize the survey weather conditions generator.
        self.w = weatherModule(self.startdate.datetime, seed)

        # Define the moonsoon season as a tuple (month, day, hour, min, sec).
        # There is no observing on monsoon_start and observations resume
        # on monsoon_stop.
        self.monsoon_start = (7, 13) + (local_noon.hour, 0, 0)
        self.monsoon_stop = (8, 27) + (local_noon.hour, 0, 0)

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
            tilesObserved.meta['MJDBEGIN'] = self.startdate.mjd
            start_val = 0
        self.tilesObserved = tilesObserved
        self.ocnt = obsCount(start_val)

        self.survey_done = False
        self.day = self.startdate
        self.iday = 0
        self.tiles_todo = self.sp.numtiles
        self.log.info(
            'Simulator initialized with {0} tiles remaining to observe.'
            .format(self.tiles_todo))


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
        assert self.day >= self.startdate and self.day < self.enddate
        self.log.info('Simulating {0}'.format(self.day.datetime.date()))

        # Lookup today's ephemerides.
        day_stats = self.surveycal[self.iday]
        sunset = astropy.time.Time(day_stats['MJDsunset'], format='mjd')
        assert sunset > self.day and sunset - self.day < 1 * u.day

        # Check if we are in the moonsoon period.
        year = self.day.datetime.year
        monsoon_start = datetime.datetime(year, *self.monsoon_start)
        monsoon_stop = datetime.datetime(year, *self.monsoon_stop)
        if (self.day.datetime < monsoon_start or
            self.day.datetime >= monsoon_stop):

            # Check if we are in the full-moon engineering period.
            if day_stats['MoonFrac'] < 0.85:

                # Simulate a normal observing night.
                ntodate = len(self.tilesObserved)
                self.w.resetDome(self.day.datetime)
                obsplan = self.sp.afternoonPlan(day_stats, self.tilesObserved)
                self.tilesObserved = nightOps(
                    day_stats, obsplan, self.w, self.ocnt, self.tilesObserved,
                    use_jpl=self.use_jpl)
                ntiles_tonight = len(self.tilesObserved)-ntodate
                self.tiles_todo -= ntiles_tonight
                self.log.info('Observed {0} tiles tonight, {1} remaining.'
                              .format(ntiles_tonight, self.tiles_todo))
                if (self.sp.numtiles - len(self.tilesObserved)) == 0:
                    self.survey_done = True
            else:
                self.log.info(
                    'No observing around full moon ({0:.1f}% illuminated).'
                    .format(100 * day_stats['MoonFrac']))
        else:
            self.log.info('No observing during monsoon {0} to {1}'
                     .format(monsoon_start.date(), monsoon_stop.date()))

        self.day += 1 * u.day
        self.iday += 1
        if self.day == self.enddate:
            self.survey_done = True

        return not self.survey_done
