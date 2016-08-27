import numpy as np
import astropy.io.fits as pyfits
#from surveysim.utils import angsep
from surveysim.exposurecalc import airMassCalculator
from surveysim.avoidobject import avoidObject, moonLoc
from surveysim.utils import mjd2lst
from datetime import datetime
from astropy.time import Time

# This is a VERY simplfied version of the next field seclector.
# The only thing it does is return the first target on the list
# that is after the current time.

MAX_AIRMASS = 40.0 # doesn't do anything...

def nextFieldSelector(obsplan, mjd, conditions, tilesObserved):

    hdulist = pyfits.open(obsplan)
    tiledata = hdulist[1].data
    moonRA = hdulist[0].header['MOONRA  ']
    moonDEC = hdulist[0].header['MOONDEC ']
    tileID = tiledata['TILEID']
    tmin = tiledata['LSTMIN']
    tmax = tiledata['LSTMAX']
    explen = tiledata['MAXEXPLEN']/240.0
    ra = tiledata['RA']
    dec = tiledata['DEC']

    lst = mjd2lst(mjd)
    dt = Time(mjd, format='mjd')
    found = False
    for i in range(len(tileID)):
        t1 = tmin[i]
        t2 = tmax[i] - explen[i]

        if t1 < t2:
            if (lst > t1 and lst < t2
                and avoidObject(dt.datetime, ra[i], dec[i])
                and airMassCalculator(ra[i], dec[i], lst) < MAX_AIRMASS):
                if ( (len(tilesObserved) > 0 and tileID[i] not in tilesObserved['TILEID']) or len(tilesObserved) == 0 ):
                    found = True
                    break
        else:
            if ( ((lst > t1 and t1 <=360.0) or (lst >= 0.0 and lst < t2))
                 and avoidObject(dt.datetime, ra[i], dec[i])
                 and airMassCalculator(ra[i], dec[i], lst) < MAX_AIRMASS):
                if ( (len(tilesObserved) > 0 and tileID[i] not in tilesObserved['TILEID']) or len(tilesObserved) == 0 ):
                    found = True
                    break

    if found == True:
        tileID = tiledata['TILEID'][i]
        RA = ra[i]
        DEC = dec[i]
        program = tiledata['PROGRAM'][i]
        Ebmv = tiledata['EBV_MED'][i]
        maxLen = tiledata['MAXEXPLEN'][i]
        moonFrac = hdulist[0].header['MOONFRAC']
        DESsn2 = 100.0 # Some made-up number -> has to be the same as the reference in exposurecalc.py
        status = tiledata['STATUS'][i]
        exposure = -1.0 # Updated after observation
        obsSN2 = -1.0   # Idem
        moondist = moonLoc(dt.datetime, RA, DEC)
        target = {'tileID' : tileID, 'RA' : RA, 'DEC' : DEC, 'Program': program, 'Ebmv' : Ebmv, 'maxLen': maxLen,
                  'MoonFrac': moonFrac, 'MoonDist': moondist, 'DESsn2': DESsn2, 'Status': status, 'Exposure': exposure, 'obsSN2': obsSN2}
    else:
        target = None
    return target

        
