import numpy as np
from astropy.time import Time

class weatherModule:

#    Simulates and contains information for weather dependent variables.

    def __init__(self, datetime, seed=None):
        # Initiate random number generators with seed
        self.rn = np.random.RandomState(seed)
        self.openDome = self.simDome(datetime) # This assumes new object instantiated each night
        self.datetime = datetime

    def simDome(self, datetime):
        # returns true or false depending on whether the dome is open
        threshold = [35.24, 44.14, 27.68, 26.73, 14.22, 15.78,
                     55.92, 48.75, 29.45, 24.44, 24.86, 34.74]
        month = datetime.datetime.month
        x = self.rn.uniform()
        answer = False
        if x > 0.01*threshold[month-1]:
            answer = True
        return answer

    def resetDome(self, datetime):
        self.openDome = self.simDome(datetime)

    def simSeeing(self, datetime):
        # Returns current seeing in arcseconds
        seeing = 1.1 + self.rn.normal(0.0, 0.1)
        if seeing < 0.65:
            seeing = 0.65
        return seeing

    def simTransparency(self, datetime):
        # Returns current (linear) transparency
        # Dummy for now
        transparency = 1.0 - self.rn.lognormal(-1.25, 0.5)
        if transparency < 0.0 : transparency = 0.0
        return transparency

    def simClouds(self, datetime):
        cloudCover = self.rn.normal(0.2, 0.1)
        if cloudCover < 0.0 : cloudCover = 0.0
        if cloudCover > 1.0 : cloudCover = 1.0
        return cloudCover

    # Should be used at the begining of the night.
    def getValues(self, time=None):
        if time == None:
            time = self.datetime.mjd
        seeing = self.simSeeing(time)
        transparency = self.simTransparency(time)
        clouds = self.simClouds(time)
        weatherNow = {'Seeing': seeing, 'Transparency': transparency,
                      'OpenDome': self.openDome, 'Clouds':clouds}
        return weatherNow

    # Small variations around current values.
    def updateValues(self, conditions, time):
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
