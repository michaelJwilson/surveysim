from __future__ import print_function, division, absolute_import

from datetime import datetime

import numpy as np

import astropy.time
import astropy.table
import astropy.units as u

import desiutil.log

import desimodel.seeing

import desisurvey.config
import desisurvey.utils


# Percentage probability of the dome being closed due to weather
# during each calendar month.
dome_closed_probability = [
    35.24, 44.14, 27.68, 26.73, 14.22, 15.78,
    55.92, 48.75, 29.45, 24.44, 24.86, 34.74]


class Weather(object):
    """Simulate weather conditions affecting observations.

    Parameters
    ----------
    start_date : datetime.date or None
        Survey starts on the evening of this date. Use the ``first_day``
        config parameter if None (the default).
    stop_date : datetime.date or None
        Survey stops on the morning of this date. Use the ``last_day``
        config parameter if None (the default).
    time_step : astropy.units.Quantity
        Time step for calculating updates. Must evenly divide 24 hours.
    seed : int or None
        Random number seed to use for simulating random processes.
    """
    def __init__(self, start_date=None, stop_date=None, time_step=5 * u.min,
                 seed=123):
        self.log = desiutil.log.get_logger()
        config = desisurvey.config.Configuration()
        self.gen = np.random.RandomState(seed)

        # Use our config to set any unspecified dates.
        if start_date is None:
            start_date = config.first_day()
        if stop_date is None:
            stop_date = config.last_day()
        num_nights = (stop_date - start_date).days
        if num_nights <= 0:
            raise ValueError('Expected start_date < stop_date.')

        # Check that the time step evenly divides 24 hours.
        steps_per_day = int(round((1 * u.day / time_step).to(1).value))
        if not np.allclose((steps_per_day * time_step).to(u.day).value, 1.):
            raise ValueError(
                'Requested time_step does not evenly divide 24 hours: {0}.'
                .format(time_step))

        # Calculate the number of times where we will tabulate the weather.
        num_rows = num_nights * steps_per_day
        self._table = astropy.table.Table()

        # Initialize column of MJD timestamps.
        t0 = desisurvey.utils.local_noon_on_date(start_date)
        times = t0 + np.arange(num_rows) * time_step
        self._table['mjd'] = times.mjd.astype(np.float32)

        # Decide whether the dome is opened on each night.
        # We currently assume this is fixed for a whole night, but
        # tabulate the status at each time so that this could be
        # updated in future to simulate partial-night weather outages.
        self._table['dome'] = np.ones(num_rows, bool)
        for i in range(num_nights):
            ij = i * steps_per_day
            month = times[ij].datetime.month
            if 100 * self.gen.uniform() < dome_closed_probability[month - 1]:
                self._table['dome'][ij:ij + steps_per_day] = False

        # Sample random seeing values.
        self._table['seeing'] = desimodel.seeing.sample(
            num_rows, 24 * 3600. / steps_per_day)

        # Sample transparency.
        # ...

        self.start_date = start_date
        self.stop_date = stop_date
        self.num_nights = num_nights
        self.steps_per_day = steps_per_day


class weatherModule:
    """
    Simulates and contains information for weather dependent variables.  It is instantiated
    at the beginning and, if applicable, at every restart of the survey simulation.
    """

    def __init__(self, dt, seed=None):
        """
        Initiate random number generators with seed.

        Args:
            dt: datetime object containing the start of the simulation
            seed: integer
        """
        self.rn = np.random.RandomState(seed)
        self.openDome = self.simDome(dt)
        self.dt = dt
        self.t = astropy.time.Time(dt)

    def simDome(self, dt):
        """
        returns true or false depending on whether the dome is open

        Args:
            dt: datetime object containing the current time

        Returns:
            bool: True if dome is open.
        """
        month = dt.month
        x = self.rn.uniform()
        answer = False
        if x > 0.01*dome_closed_probability[month-1]:
            answer = True
        return answer

    def resetDome(self, dt):
        """
        Checks of the dome can open or not

        Args:
            dt: datetime object containing the current time
        """
        self.openDome = self.simDome(dt)

    def simSeeing(self, dt):
        """
        Draws the seeing from lognormal distribution.
        Approximates the results from Dey & Valdes.

        Args:
            dt: datetime object containing the current time

        Returns:
            seeing: float (arcseconds)
        """
        seeing = self.rn.lognormal(0.0, 0.25)
        if seeing < 0.65:
            seeing = 0.65
        return seeing

    def simTransparency(self, dt):
        """
        Computes current (linear) transparency, drawn from lognormal distribution.

        Args:
            dt: datetime object containing the current time

        Returns:
            transparency: float
        """
        transparency = self.rn.lognormal(0.11111111, 0.3333333)
        if transparency < 0.0 : transparency = 0.0
        if transparency > 1.0 : transparency = 1.0
        return transparency

    def simClouds(self, dt):
        """
        Draws cloud cover from a normal distribution

        Args:
            dt: datetime object containing the current time

        Returns:
            cloudCover: float, between 0 and 1
        """
        cloudCover = self.rn.normal(0.2, 0.1)
        if cloudCover < 0.0 : cloudCover = 0.0
        if cloudCover > 1.0 : cloudCover = 1.0
        return cloudCover

    def getValues(self, time=None):
        """
        Gets all weather values; should be called at the begining of the night.

        Args:
            dt: datetime object containing the current time

        Returns:
            dictionnary containing the follwing keys:
            'Seeing', 'Transparency', 'OpenDome', 'Clouds'
        """
        if time == None:
            time = self.dt.mjd
        seeing = self.simSeeing(time)
        transparency = self.simTransparency(time)
        clouds = self.simClouds(time)
        weatherNow = {'Seeing': seeing, 'Transparency': transparency,
                      'OpenDome': self.openDome, 'Clouds':clouds}
        return weatherNow

    def updateValues(self, conditions, time):
        """
        Computes small variations around current weather values.
        Leaves dome as it is.

        Args:
            dt: datetime object containing the current time

        Returns:
            dictionnary containing the follwing keys:
            'Seeing', 'Transparency', 'OpenDome', 'Clouds'
        """
        seeing = conditions['Seeing'] + self.rn.normal(0.0, 0.0167)
        if seeing < 0.65:
            seeing = 0.65
        transparency = conditions['Transparency'] + self.rn.normal(0.0, 0.02)
        if transparency < 0.0:
            transparency = 0.0
        if transparency > 1.0:
            transparency = 1.0
        clouds = conditions['Clouds'] + self.rn.normal(0.0, 0.05)
        if clouds < 0.0:
            clouds = 0.0
        if clouds > 1.0:
            clouds = 1.0
        weatherNow = {'Seeing': seeing, 'Transparency': transparency,
                      'OpenDome': self.openDome, 'Clouds':clouds}
        return weatherNow
