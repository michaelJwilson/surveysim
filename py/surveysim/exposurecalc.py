import numpy as np
from astropy.time import Time
from surveysim.weather import weatherModule
from surveysim.utils import radec2altaz

def expTimeEstimator(weatherNow, amass, program, ebmv, sn2, moonFrac):
    # Estimates expusure length given current conditions.

    seeing_ref = 1.1 # Seeing value to which actual seeing is normalised
    exp_ref_dark = 1000.0   # Reference exposure time in seconds
    exp_ref_bright = 300.0  # Idem but for bright time programme
    exp_ref_grey = exp_ref_dark
    sn2_nom = 100.0 # Nominal sign-to-noise: again, made-up number

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
    f_seeing =  (a+b*seeing+c*seeing*seeing) / (a-0.25*b*b/c)
    if weatherNow['Transparency'] > 0.0:
        f_transparency = 1.0 / weatherNow['Transparency']
    else:
        f_transparency = 1.0e9
    """
    Ag=3.303*ebv[i]
    Ai=1.698*ebv[i]
    i_increase[i]=(10^(Ai/2.5))^2
    g_increase[i]=(10^(Ag/2.5))^2
    """
    Ag = 3.303*ebmv # Use g-band
    f_ebmv = np.power(10.0,Ag/2.5)
    f_am = np.power(amass,1.25)
    """
    if moonFrac < 1.0:
        f_moon = 1.0 / (1.0 - moonFrac/100.0)
    else:
        f_moon = 30.0
    """
    f_moon = 1.0 # Temporary until real values are in the code
    #print (f_am, f_seeing, f_transparency, f_ebmv, f_moon)
    f = f_am * f_seeing * f_transparency * f_ebmv * f_moon
    if f >= 0.0:
        value = exp_ref * f
    else:
        value = exp_ref
    return value

def airMassCalculator(ra, dec, lst): # Valid for small to moderate angles.
    Alt, Az = radec2altaz(ra, dec, lst)
    # Rosenberg (1966) formula
    cosZ = np.cos(np.radians(90.0-Alt))
    if Alt >= 0.0:
        amass = 1.0/(cosZ + 0.025*np.exp(-11.0*cosZ))
        if amass <= 0.0:
            print ('ERROR: negative airmass (', amass, '); LST, RA, DEC = ', lst, ra, dec)
            print ('Alt, Az = ', Alt, Az)
    else:
        amass = 1.0e99

    return amass
