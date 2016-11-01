import ephem
from datetime import datetime
from numpy import pi as PI
from surveysim.kpno import mayall

MIN_VENUS_SEP = 2.0 * PI / 180.0
MIN_MARS_SEP = 2.0 * PI / 180.0
MIN_JUPITER_SEP = 2.0 * PI / 180.0
MIN_SATURN_SEP = 2.0 * PI / 180.0
MIN_NEPTUNE_SEP = 2.0 * PI / 180.0
MIN_URANUS_SEP = 2.0 * PI / 180.0
MIN_CERES_SEP = 2.0 * PI / 180.0

def avoidObject(datetime, ra0, dec0):
    """
    Returns True if all the objects on the list are far enough away from
    the coordinates (assumed to be apparent or observed, not mean).  The datetime
    object should have timezone info included.  The inputs are in decimal degrees.

    The current list has: Venus, Mars, Jupiter, Saturn, Neptune, Uranus;
    the Moon is treated separately.
    """

    ra = PI * ra0/180.0
    dec = PI * dec0/180.0
    
    dt = ephem.Date(datetime)
    gatech = ephem.Observer()
    gatech.lon, gatech.lat = mayall.west_lon_deg, mayall.lat_deg
    gatech.date = dt
    gatech.epoch = dt

    venus = ephem.Venus()
    venus.compute(gatech)
    if ephem.separation(venus, (ra, dec)) < MIN_VENUS_SEP:
        return False
    mars = ephem.Mars()
    mars.compute(gatech)
    if ephem.separation(mars, (ra, dec)) < MIN_MARS_SEP:
        return False
    #ceres = ephem.Ceres()
    #ceres.compute(gatech)
    #if ephem.separation(ceres, (ra, dec)) < MIN_CERES_SEP:
    #    return False
    jupiter = ephem.Jupiter()
    jupiter.compute(gatech)
    if ephem.separation(jupiter, (ra, dec)) < MIN_JUPITER_SEP:
        return False
    saturn = ephem.Saturn()
    saturn.compute(gatech)
    if ephem.separation(saturn, (ra, dec)) < MIN_SATURN_SEP:
        return False
    neptune = ephem.Neptune()
    neptune.compute(gatech)
    if ephem.separation(neptune, (ra, dec)) < MIN_NEPTUNE_SEP:
        return False
    uranus = ephem.Uranus()
    uranus.compute(gatech)
    if ephem.separation(uranus, (ra, dec)) < MIN_URANUS_SEP:
        return False

    # If still here, return True
    return True

def moonLoc (datetime, ra0, dec0):
    """
    Returns the distance to the Moon if RA and DEC as well as alt, az.
    Input and outputs are 
    """

    dt = ephem.Date(datetime)
    gatech = ephem.Observer()
    gatech.lon, gatech.lat = mayall.west_lon_deg, mayall.lat_deg
    gatech.date = dt
    gatech.epoch = dt

    moon = ephem.Moon()
    moon.compute(gatech)
    ra = PI * ra0/180.0
    dec = PI * dec0/180.0
    moondist = ephem.separation(moon, (ra, dec))
    moondist*180.0/PI
    
    return moondist*180.0/PI, (moon.alt)*180.0/PI, (moon.az)*180.0/PI
