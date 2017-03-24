import numpy as np
from astropy.time import Time
from datetime import datetime

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
        print(dt)
        self.rn = np.random.RandomState(seed)
        self.openDome = self.simDome(dt)
        self.dt = dt
        self.t = Time(dt)

    def simDome(self, dt):
        """
        returns true or false depending on whether the dome is open

        Args:
            dt: datetime object containing the current time

        Returns:
            bool: True if dome is open.
        """
        threshold = [35.24, 44.14, 27.68, 26.73, 14.22, 15.78,
                     55.92, 48.75, 29.45, 24.44, 24.86, 34.74]
        month = dt.month
        x = self.rn.uniform()
        answer = False
        if x > 0.01*threshold[month-1]:
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
