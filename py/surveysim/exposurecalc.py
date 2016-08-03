import numpy as np
from astropy.time import Time
from surveysim.weather import weatherModule
from surveysim.utils import radec2altaz

def expTimeEstimator(weatherNow, airmass, program, ebmv, sn2, moonFrac):
# Estimates expusure length given current conditions.

    seeing_ref = 1.1 # Seeing value to which actual seeing is normalised
    exp_ref_dark = 1000.0   # Reference exposure time in seconds
    exp_ref_bright = 300.0  # Idem but for bright time programme
    exp_ref_grey = 650.0    # Made up number: just took the average
    sn2_nom = 100.0 # Nominal sign-to-noise

    if program == "DARK":
        exp_ref = exp_ref_dark
    elif program == "BRIGHT":
        exp_ref = exp_ref_bright
    elif program == "GRAY":
        exp_ref = exp_ref_grey
    else:
        exp_ref = 0.0 # Replace with throwing an exception
    seeing = weatherNow['Seeing']
    a = 4.6
    b = -1.55
    c = 1.15
    f_seeing = (a-0.25*b*b/c) / (a+b*seeing+c*seeing*seeing)
    #print (weatherNow['Transparency'])
    if weatherNow['Transparency'] > 0.0:
        f_transparency = 1.0 / weatherNow['Transparency']
    else:
        f_transparency = 1.0e9
    f_ebmv = np.power(10.0,ebmv/2.5)
    """
    if moonFrac < 1.0:
        f_moon = 1.0 / (1.0 - moonFrac/100.0)
    else:
        f_moon = 30.0
    """
    f_moon = 1.0 # Temporary until real values are in the code
    f = f_seeing * f_transparency * f_ebmv * f_moon
    if f >= 0.0:
        value = exp_ref * f * (sn2 / sn2_nom)
    else:
        value = exp_ref
    return value

def airMassCalculator(ra, dec, lst): # Valid for small to moderate angles.
    Alt, Az = radec2altaz(ra, dec, lst)
    # Rosenberg (1966) formula
    cosZ = np.cos(np.radians(90.0-Alt))
    amass = 1.0/(cosZ + 0.025*np.exp(-11.0*cosZ))
    if amass <= 0.0:
        print ('ERROR: negative airmass (', amass, '); LST, RA, DEC = ', lst, ra, dec,'!')
    return amass
