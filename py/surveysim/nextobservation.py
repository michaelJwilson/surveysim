import numpy as np
import astropy.io.fits as pyfits

# This is a VERY simplfied version of the next field seclector.
# The only thing it does is return the first target on the list
# that is after the current time.

MIN_MOON_SEP = 5.0   # In degrees

def nextFieldSelector(obsplan, lst, conditions, tilesObserved):

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

    found = False
    i = 0
    for t1 in np.nditer(tmin):
        t2 = tmax[i] - explen[i]
        if (lst > t1 and lst < t2
            and abs(ra[i] - moonRA) > MIN_MOON_SEP and abs(dec[i] - moonDEC) > MIN_MOON_SEP
            and tileID[i] not in tilesObserved):
            found = True
            break
        i = i+1

    if found == True:
        tileID = tiledata['TILEID'][i]
        RA = tiledata['RA'][i]
        DEC = tiledata['DEC'][i]
        program = 'Dark' # Needs to be in the list
        Ebmv = tiledata['EBV_MED'][i]
        maxLen = tiledata['MAXEXPLEN'][i]
        moonFrac = hdulist[0].header['MOONFRAC']
        DESsn2 = 100.0 # Some made-up number
        status = tiledata['STATUS'][i]
        exposure = -1.0 # Updated after observation
        obsSN2 = -1.0   # Idem

        target = {'tileID' : tileID, 'RA' : RA, 'DEC' : DEC, 'Program': program, 'Ebmv' : Ebmv, 'maxLen': maxLen,
                  'MoonFrac': moonFrac, 'DESsn2': DESsn2, 'Status': status, 'Exposure': exposure, 'obsSN2': obsSN2}
    else:
        target = None
    return target

        
