import numpy as np
from astropy.time import Time
import astropy.units as u
from surveysim.weather import weatherModule
from surveysim.utils import radec2altaz

def expTimeEstimator(weatherNow, amass, program, ebmv, sn2, moonFrac, moonDist, moonAlt):
    """
    Estimates expusure length given current conditions.

    Args:
        weatherNow: dictionnary containing the following keys:
                    'Seeing', 'Transparency', 'OpenDome', 'Clouds'
        amass: float, air mass
        programm: string, 'DARK', 'BRIGHT' or 'GRAY'
        ebmv: float, E(B-V)
        sn2: float, desired (S/N)^2
        moonFrac: float, Moon illumination fraction, between 0 (new) and 1 (full).
        moonDist: float, separation angle between field center and moon in degrees.
        moonAlt: float, moon altitude angle in degrees.

    Returns:
        float, estimated exposure time
    """

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

    #Ag=3.303*ebv[i]
    #Ai=1.698*ebv[i]
    #i_increase[i]=(10^(Ai/2.5))^2
    #g_increase[i]=(10^(Ag/2.5))^2

    Ag = 3.303*ebmv # Use g-band
    f_ebmv = np.power(10.0,Ag/2.5)
    f_am = np.power(amass,1.25)

    f_moon = moonExposureTimeFactor(moonFrac, moonDist, moonAlt)
    #print (f_am, f_seeing, f_transparency, f_ebmv, f_moon)
    f = f_am * f_seeing * f_transparency * f_ebmv * f_moon
    if f >= 0.0:
        value = exp_ref * f
    else:
        value = exp_ref
    return value


def moonExposureTimeFactor(moonFrac, moonDist, moonAlt):
    """Calculate exposure time factor due to scattered moonlight.

    This factor is based on a study of SNR for ELG targets and designed to
    achieve a median SNR of 7 for a typical ELG [OII] doublet at the lower
    flux limit of 8e-17 erg/(cm2 s A), averaged over the expected ELG target
    redshift distribution 0.6 < z < 1.7.

    For details, see the jupyter notebook doc/nb/ScatteredMoon.ipynb in
    this package.

    Parameters
    ----------
    moonFrac : float
        Illuminated fraction of the moon, between 0-1.
    moonDist : float
        Separation angle between field center and moon in degrees.
    moonAlt : float
        Altitude angle of the moon above the horizon in degrees.

    Returns
    -------
    float
        Dimensionless factor that exposure time should be increased to
        account for increased sky brightness due to scattered moonlight.
        Will be 1 when the moon is below the horizon.
    """
    return 1.


def airMassCalculator(ra, dec, lst): # Valid for small to moderate angles.
    """
    Calculates airmass given position and LST.  Uses formula from
    Rosenberg (1966)

    Args:
        ra: float (degrees)
        dec: float (degrees)
        lst: float (degrees)

    Returns:
        float, air mass
    """

    Alt, Az = radec2altaz(ra, dec, lst)
    cosZ = np.cos(np.radians(90.0-Alt))
    if Alt >= 0.0:
        amass = 1.0/(cosZ + 0.025*np.exp(-11.0*cosZ))
        if amass <= 0.0:
            print ('ERROR: negative airmass (', amass, '); LST, RA, DEC = ', lst, ra, dec)
            print ('Alt, Az = ', Alt, Az)
    else:
        amass = 1.0e99

    return amass
