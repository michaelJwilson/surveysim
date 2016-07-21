import numpy as np
from astropy.time import Time
from weather import weatherModule

# Location of Kitt Peak National Observatory Mayall 4m Telescope
# Taken from SlaLib obs.c, which cites the 1981 Almanac
Lon_KPNO_deg = -1.0 * (111.0 + (35.0 + 57.61/60.0)/60.0)
Lat_KPNO_deg = 31.0 + (57.0 + 50.3/60.0)/60.0
Alt_KPNO_m   = 2120.0

def expTimeEstimator(weatherNow, airmass, program, ebmv, sn2, moonFrac):
# Estimates expusure length given current conditions.

    seeing_ref = 1.1 # Seeing value to which actual seeing is normalised
    exp_ref_dark = 1000.0   # Reference exposure time in seconds
    exp_ref_bright = 300.0  # Idem but for bright time programme
    f_ref = 0.85
    sn2_nom = 100.0 # Nominal sign-to-noise

    if program == "Dark":
        exp_ref = exp_ref_dark
    elif program == "Bright":
        exp_ref = exp_ref_bright
    else:
        exp_ref = 0.0 # Replace with throwing an exception
    f_seeing = weatherNow['Seeing'] / seeing_ref
    f_transparency = weatherNow['Transparency']
    f_airmass = 1.0/airmass
    f_ebmv = np.exp(-ebmv) # What's the correct factor?
    f_moon = 1.0 - moonFrac/100.0 # What is the correct attenuation due to the increased sky brightness?
    f = f_seeing * f_transparency * f_airmass * f_ebmv * f_moon
    if f > 0.0:
        value = exp_ref / (f/f_ref) * (sn2 / sn2_nom)
    else:
        value = exp_ref
    return value

def airMassCalculator(ra, dec, lst): # Valid for small to moderate angles.
    # RA and LST are in decimal hours, DEC in decimal degrees
    # Note that these are *observed* RA and DEC, not mean, not apparent.
    h = np.radians(lst - ra)
    d = np.radians(dec)
    phi = np.radians(Lat_KPNO_deg)

    # sinAz = np.sin(h) / (np.cos(h)*np.sin(phi) - np.tan(d)*np.cos(phi))
    sina = np.sin(phi)*np.sin(d) + np.cos(phi)*np.cos(d)*np.cos(h)

    amass = 1.0/sina
    if amass <= 0.0:
        print 'ERROR: negative airmass (', amass, ') !!!'
    return amass

