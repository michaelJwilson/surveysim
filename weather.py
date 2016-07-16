import numpy as np
from astropy.time import Time

class weatherModule:

#    Simulates and contains information for weather dependent variables.

    def __init__(self, datetime, seeds=None):
        # Initiate random number generators with seeds
        self.openDome = self.simDome(datetime) # This assumes new object instantiated each night
        self.datetime = datetime

    def simDome(self, datetime):
        # returns true or false depending on whether the dome is open
        threshold = [35.24, 44.14, 27.68, 26.73, 14.22, 15.78,
                     55.92, 48.75, 29.45, 24.44, 24.86, 34.74]
        a = str(datetime.datetime)
        month = int(a[5]) * 10 + int(a[6])
        x = np.random.random(None)
        answer = False
        if x > 0.01*threshold[month-1]:
            answer = True
        return answer

    def resetDome(self, datetime):
        self.openDome = self.simDome(datetime)

    def simSeeing(self, datetime):
        # Returns current seeing in arcseconds
        # Dummy for now
        seeing = 1.1 + np.random.normal(0.0, 0.1)
        return seeing

    def simTransparency(self, datetime):
        # Returns current (linear) transparency
        # Dummy for now
        transparency = 1.0 - np.random.lognormal(-1.25, 0.5)
        if transparency < 0.0 : transparency = 0.0
        return transparency

    def simClouds(self, datetime):
        cloudCover = np.random.normal(0.2, 0.1)
        if cloudCover < 0.0 : cloudCover = 0.0
        if cloudCover > 1.0 : cloudCover = 1.0
        return cloudCover


    def getValues(self, time=None):
        if time == None:
            time = self.datetime.mjd
        seeing = self.simSeeing(time)
        transparency = self.simTransparency(time)
        clouds = self.simClouds(time)
        weatherNow = {'Seeing': seeing, 'Transparency': transparency,
                      'OpenDome': self.openDome, 'Clouds':clouds}
        return weatherNow
