import numpy as np
from astropy.time import Time

# Location of Kitt Peak National Observatory Mayall 4m Telescope
# Taken from SlaLib obs.c, which cites the 1981 Almanac
Lon_KPNO_deg = -1.0 * (111.0 + (35.0 + 57.61/60.0)/60.0)
Lat_KPNO_deg = 31.0 + (57.0 + 50.3/60.0)/60.0
Alt_KPNO_m   = 2120.0

# Converts decimal MJD to LST in decimal degrees
def mjd2lst(mjd):
    t = Time(mjd, format = 'mjd', location=('-111.6d', '32.0d'))
    lst_tmp = t.copy()
    try:
        lst_str = str(lst_tmp.sidereal_time('apparent'))
    except IndexError:
        lst_tmp.delta_ut1_utc = -0.1225
        lst_str = str(lst_tmp.sidereal_time('apparent'))
        # 23h09m35.9586s
        # 01234567890123
    if lst_str[2] == 'h':
        lst_hr = float(lst_str[0:2])
        lst_mn = float(lst_str[3:5])
        lst_sc = float(lst_str[6:-1])
    else:
        lst_hr = float(lst_str[0:1])
        lst_mn = float(lst_str[2:4])
        lst_sc = float(lst_str[5:-1])
    lst = lst_hr + lst_mn/60.0 + lst_sc/3600.0
    lst *= 15.0 # Convert from hours to degrees
    return lst

# All quantities are in decimal degrees
# Note that these should be *observed* RA and DEC, not mean, not apparent.
def radec2altaz(ra, dec, lst):
    
    h = np.radians(lst - ra)
    if h < 0.0:
        h += 360.0
    d = np.radians(dec)
    phi = np.radians(Lat_KPNO_deg)

    sinAz = np.sin(h) / (np.cos(h)*np.sin(phi) - np.tan(d)*np.cos(phi))
    sinAlt = np.sin(phi)*np.sin(d) + np.cos(phi)*np.cos(d)*np.cos(h)

    if sinAlt > 1.0:
        sinAlt = 1.0
    if sinAlt < -1.0:
        sinAlt = -1.0
    if sinAz > 1.0:
        sinAz = 1.0
    if sinAz < -1.0:
        sinAz = -1.0

    return np.arcsin(sinAlt), np.arcsin(sinAz)

# Calculates the angular separation between two objects.
# All quantities are in decimal degrees.
def angsep(ra1, dec1, ra2, dec2):
    deltaRA = np.radians(ra1-ra2)
    DEC1 = np.radians(dec1)
    DEC2 = np.radians(dec2)
    cosDelta = np.sin(DEC1)*np.sin(DEC2) + np.cos(DEC1)*np.cos(DEC2)*np.cos(deltaRA)
    return np.degrees(np.arccos(cosDelta))
