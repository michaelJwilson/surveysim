from astropy.time import Time
from datetime import datetime
import numpy as np
import ephem

# Assumes that the input time is midday local time
def getCal(day):

    mayall = ephem.Observer()
    mayall.lat, mayall.lon = (31.9639722222)/180.0*np.pi, -111.59933611/180.0*np.pi
    mayall.date = day
    mayall.pressure = 0.0      # Disabling pyephem's refraction calculations, just use
    mayall.horizon = '-0:34'   # the value that the USNO uses.
    day0 = Time(datetime(1900,1,1,12,0,0)) # This throws a warning because of the early year, but it is harmless.
    # Sun
    MJDsunset = float( mayall.next_setting(ephem.Sun()) + day0.mjd )
    MJDsunrise = float( mayall.next_rising(ephem.Sun()) + day0.mjd )
    # 13 deg twilight, adequate (?) for BGS sample.
    mayall.horizon = '-13'
    MJDe13twi = float( mayall.next_setting(ephem.Sun(), use_center=True) + day0.mjd )
    MJDm13twi = float( mayall.next_rising(ephem.Sun(), use_center=True) + day0.mjd )
    # 15 deg twilight, start of dark time if the moon is down.
    mayall.horizon = '-15'
    MJDetwi = float( mayall.next_setting(ephem.Sun(), use_center=True) + day0.mjd )
    MJDmtwi = float( mayall.next_rising(ephem.Sun(), use_center=True) + day0.mjd )   
    # Moon.
    mayall.horizon = '-0:34'
    MJDmoonrise = float( mayall.next_rising(ephem.Moon()) + day0.mjd )
    MJDmoonset = float( mayall.next_setting(ephem.Moon()) + day0.mjd )
    if (MJDmoonrise > MJDsunrise):
       MJDmoonrise = float( mayall.previous_rising(ephem.Moon()) + day0.mjd )
    m0 = ephem.Moon()
    m0.compute(day)
    MoonFrac = float( m0.moon_phase )

    # Get the night's directory name right away.
    if day.month >= 10:
        monthStr = str(day.month)
    else:
        monthStr = '0' + str(day.month)
    if day.day >= 10:
        dayStr = str(day.day)
    else:
        dayStr = '0' + str(day.day)
    dirName = str(day.year) + monthStr + dayStr
    day_stats = {'MJDsunset': MJDsunset,
                 'MJDsunrise': MJDsunrise,
                 'MJDetwi': MJDetwi,
                 'MJDmtwi': MJDmtwi,
                 'MJDe13twi': MJDe13twi,
                 'MJDm13twi': MJDm13twi,
                 'MJDmoonrise': MJDmoonrise,
                 'MJDmoonset': MJDmoonset,
                 'MoonFrac' : MoonFrac,
                 'dirName': dirName}
    return day_stats

