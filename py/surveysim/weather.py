import numpy as np
from astropy.time import Time
from datetime import datetime

class weatherModule:

    # Simulates and contains information for weather dependent variables.

    def __init__(self, dt, seed=None):
        # Initiate random number generators with seed
        print(dt)
        self.rn = np.random.RandomState(seed)
        self.openDome = self.simDome(dt) # This assumes new object instantiated each night
        self.dt = dt
        self.t = Time(dt)

    def simDome(self, dt):
        # returns true or false depending on whether the dome is open
        threshold = [35.24, 44.14, 27.68, 26.73, 14.22, 15.78,
                     55.92, 48.75, 29.45, 24.44, 24.86, 34.74]
        month = dt.month
        x = self.rn.uniform()
        answer = False
        if x > 0.01*threshold[month-1]:
            answer = True
        return answer

    def resetDome(self, dt):
        self.openDome = self.simDome(dt)

    def simSeeing(self, dt):
        # Returns current seeing in arcseconds
        seeing = self.rn.lognormal(0.0, 0.25)
        if seeing < 0.5:
            seeing = 0.5
        return seeing

    def simTransparency(self, dt):
        # Returns current (linear) transparency
        transparency = self.rn.lognormal(0.11111111, 0.3333333)
        if transparency < 0.0 : transparency = 0.0
        if transparency > 1.0 : transparency = 1.0
        return transparency

    def simClouds(self, dt):
        cloudCover = self.rn.normal(0.2, 0.1)
        if cloudCover < 0.0 : cloudCover = 0.0
        if cloudCover > 1.0 : cloudCover = 1.0
        return cloudCover

    # Should be used at the begining of the night.
    def getValues(self, time=None):
        if time == None:
            time = self.dt.mjd
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
