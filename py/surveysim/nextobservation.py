import numpy as np
import astropy.io.fits as pyfits
#from surveysim.utils import angsep
from surveysim.exposurecalc import airMassCalculator
from surveysim.avoidobject import avoidObject, moonLoc
from surveysim.utils import mjd2lst
from datetime import datetime
from astropy.time import Time

MAX_AIRMASS = 10.0 #3.0 This new bound effectively does nothing.
MIN_MOON_SEP = 90.0

def nextFieldSelector(obsplan, mjd, conditions, tilesObserved):
    """
    Returns the first tile for which the current time falls inside
    its assigned LST window and is far enough from the Moon and
    planets.

    Args:
        obsplan: string, FITS file containing the afternoon plan
        mjd: float, current time
        conditions: dictionnary containing the weather info
        tilesObserved: list containing the tileID of all completed tiles

    Returns:
        target: dictionnary containing the following keys:
                'tileID', 'RA', 'DEC', 'Program', 'Ebmv', 'maxLen',
                'MoonFrac', 'MoonDist', 'MoonAlt', 'DESsn2', 'Status',
                'Exposure', 'obsSN2'
    """

    hdulist = pyfits.open(obsplan)
    tiledata = hdulist[1].data
    moonfrac = hdulist[0].header['MOONFRAC']
    tileID = tiledata['TILEID']
    tmin = tiledata['LSTMIN']
    tmax = tiledata['LSTMAX']
    explen = tiledata['MAXEXPLEN']/240.0
    ra = tiledata['RA']
    dec = tiledata['DEC']
    program = tiledata['PROGRAM']

    lst = mjd2lst(mjd)
    dt = Time(mjd, format='mjd')
    found = False
    for i in range(len(tileID)):
        t1 = tmin[i]
        t2 = tmax[i] - explen[i]

        if ( ((t1 <= t2) and (lst > t1 and lst < t2)) or ( (t2 < t1) and ((lst > t1 and t1 <=360.0) or (lst >= 0.0 and lst < t2))) ):
            if (avoidObject(dt.datetime, ra[i], dec[i]) and airMassCalculator(ra[i], dec[i], lst) < MAX_AIRMASS):
                moondist, moonalt, moonaz = moonLoc(dt.datetime, ra[i], dec[i])
                if ( (len(tilesObserved) > 0 and tileID[i] not in tilesObserved['TILEID']) or len(tilesObserved) == 0 ):
                    if ( (moonalt < 0.0 and program[i] == 'DARK') or
                         ((moonfrac < 0.2 or (moonalt*moonfrac < 12.0)) and moondist > MIN_MOON_SEP and program[i] == 'GRAY') or
                         (program[i] == 'BRIGHT') ):
                        found = True
                        break

    if found == True:
        tileID = tiledata['TILEID'][i]
        RA = ra[i]
        DEC = dec[i]
        Ebmv = tiledata['EBV_MED'][i]
        maxLen = tiledata['MAXEXPLEN'][i]
        DESsn2 = 100.0 # Some made-up number -> has to be the same as the reference in exposurecalc.py
        status = tiledata['STATUS'][i]
        exposure = -1.0 # Updated after observation
        obsSN2 = -1.0   # Idem
        target = {'tileID' : tileID, 'RA' : RA, 'DEC' : DEC, 'Program': program[i], 'Ebmv' : Ebmv, 'maxLen': maxLen,
                  'MoonFrac': moonfrac, 'MoonDist': moondist, 'MoonAlt': moonalt, 'DESsn2': DESsn2, 'Status': status, 'Exposure': exposure, 'obsSN2': obsSN2}
    else:
        target = None
    return target

