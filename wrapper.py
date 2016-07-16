#! /usr/bin/python

from astropy.time import Time
from datetime import datetime
from weather import weatherModule
from calendar import obsCalendar
from exposurecalc import expTimeEstimator
from exposurecalc import airMassCalculator
import astropy.io.fits as pyfits
from afternoonplan import surveyPlan
from nextobservation import nextFieldSelector
from observefield import observeField

def nightOps(day_stats, obsplan, w):

    nightOver = False
    mjd = day_stats['MJDsunset']
    tilesObserved = []
    tileIDdone = []

    conditions = w.getValues(mjd)
    print conditions
    if conditions['OpenDome'] == True:

        while nightOver == False:
            conditions = w.getValues(mjd)

            # t = Time(mjd, format = 'mjd', location=('-111.6d', '32.0d'))
            # lst = t.sidereal_time('apparent')
            lst = (mjd - float(int(mjd)) - 7.0/24.0) * 360.0 # temporary hack
            # print 'LST = ', lst
            target = nextFieldSelector(obsplan, lst, conditions, tileIDdone)
            if target != None:
                # Compute mean to apparent to observed ra and dec???
                airmass = airMassCalculator(target['RA'], target['DEC'], lst)
                exposure = expTimeEstimator(conditions, airmass, target['Program'], target['Ebmv'], target['DESsn2'], day_stats['MoonFrac'])
                # print 'Estimated exposure = ', exposure, 'Maximum allowed exposure = ', target['maxLen']
                if exposure < target['maxLen']:
                    status, real_exposure, real_sn2 = observeField(target, exposure)
                    target['Status'] = status
                    target['exposure'] = real_exposure
                    target['obsSN2'] = real_sn2
                    mjd += real_exposure/86400.0
                    tilesObserved.append(target)
                    tileIDdone.append(target['tileID']) # This is for an easy check in the nextFieldSelector.
                    # Output headers, but no data.
                    # In the future: GFAs (i, x, y + metadata for i=id, time, postagestampid) and fiber positions.
                    prihdr = pyfits.Header()
                    prihdr['TILEID  '] = target['tileID']
                    prihdr['RA      '] = target['RA']
                    prihdr['DEC     '] = target['DEC']
                    prihdr['PROGRAM '] = target['Program']
                    prihdr['EBMV    '] = target['Ebmv']
                    prihdr['MAXLEN  '] = target['maxLen']
                    prihdr['MOONFRAC'] = target['MoonFrac']
                    prihdr['DESSN2  '] = target['DESsn2']
                    prihdr['STATUS  '] = target['Status']
                    prihdr['EXPLEN  '] = target['Exposure']
                    prihdr['OBSSN2  '] = target['obsSN2']
                    t = Time(mjd, format = 'mjd')
                    prihdr['DATE-OBS'] = t.fits
                    filename = 'tileID' + str(target['tileID']) + '_' + str(t.isot) + '.fits'
                    prihdu = pyfits.PrimaryHDU(header=prihdr)
                    prihdu.writeto(filename, clobber=True)
                else:
                    # Choose another target?
                    # Observe longer split into modulo(max_len)
                    mjd += 0.5/24.0
            else:
                mjd += 0.5/24.0
            # Check time
            if mjd > day_stats['MJDsunrise']:
                nightOver = True
        
    return tilesObserved

def surveySim(startday, startmonth, startyear):

    sp = surveyPlan()
    startday = Time(datetime(startyear, startmonth, startday, 12, 0, 0))
    mjd_start = startday.mjd
    w = weatherModule(startday)

    cal = obsCalendar("calendar.2016")
    iday = 0
    for day in cal:
        t = Time(day['MJDsunset'], format = 'mjd')
        w.resetDome(t)
        if day['MJDsunset'] < mjd_start:
            continue
        if iday == 0:
            tiles_todo, obsplan = sp.afternoonPlan(day)
        else:
            tiles_todo, obsplan = sp.afternoonPlan(day, tiles_observed)
        tiles_observed = nightOps(day, obsplan, w)
        t = Time(day['MJDsunset'], format = 'mjd')
        print 'On the night starting ', t.iso, ', we observed ', len(tiles_observed), ' tiles.'
        if tiles_todo == 0:
            break
        iday = iday+1

