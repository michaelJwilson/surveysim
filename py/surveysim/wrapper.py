#! /usr/bin/python

import numpy as np
from astropy.time import Time
from datetime import datetime
from surveysim.weather import weatherModule
from surveysim.calendar import obsCalendar
from surveysim.exposurecalc import expTimeEstimator
from surveysim.exposurecalc import airMassCalculator
import astropy.io.fits as pyfits
from surveysim.afternoonplan import surveyPlan
from surveysim.nextobservation import nextFieldSelector
from surveysim.observefield import observeField
from astropy.table import Table
import os
from surveysim.utils import mjd2lst

class obsCount:

    def __init__(self):
        self.obsNumber = 0
    
    def update(self):
        self.obsNumber += 1
        obsNumber = self.obsNumber
        if obsNumber >= 10000000:
            partFileName = str(obsNumber)
        elif obsNumber < 10000000 and obsNumber >= 1000000:
            partFileName = '0' + str(obsNumber)
        elif obsNumber < 1000000 and obsNumber >= 100000:
            partFileName = '00' + str(obsNumber)
        elif obsNumber < 100000 and obsNumber >= 10000:
            partFileName = '000' + str(obsNumber)
        elif obsNumber < 10000 and obsNumber >= 1000:
            partFileName = '0000' + str(obsNumber)
        elif obsNumber < 1000 and obsNumber >= 100:
            partFileName = '00000' + str(obsNumber)
        elif obsNumber < 100 and obsNumber >= 10:
            partFileName = '000000' + str(obsNumber)
        else:
            partFileName = '0000000' + str(obsNumber)
        return partFileName

def nightOps(day_stats, obsplan, w, ocnt, tilesObserved, tableOutput=True):

    nightOver = False
    mjd = day_stats['MJDsunset']

    if tableOutput:
        obsList = []
    else:
        os.mkdir(day_stats['dirName'])
    
    conditions = w.getValues(mjd)
    print ('Conditions at the beginning of the night: ')
    print (conditions)
    if conditions['OpenDome'] == True:

        while nightOver == False:
            conditions = w.updateValues(conditions, mjd)

            lst = mjd2lst(mjd)
            target = nextFieldSelector(obsplan, lst, conditions, tilesObserved)
            if target != None:
                # Compute mean to apparent to observed ra and dec???
                airmass = airMassCalculator(target['RA'], target['DEC'], lst)
                exposure = expTimeEstimator(conditions, airmass, target['Program'], target['Ebmv'], target['DESsn2'], day_stats['MoonFrac'])
                #print ('Estimated exposure = ', exposure, 'Maximum allowed exposure for tileID', target['tileID'], ' = ', target['maxLen'])
                if exposure < target['maxLen']:
                    status, real_exposure, real_sn2 = observeField(target, exposure)
                    target['Status'] = status
                    target['Exposure'] = real_exposure
                    target['obsSN2'] = real_sn2
                    mjd += real_exposure/86400.0
                    exposureAttempted = True
                    tilesObserved.add_row([target['tileID'], status])
                    if tableOutput:
                        t = Time(mjd, format = 'mjd')
                        tbase = str(t.isot)
                        obsList.append((target['tileID'],  target['RA'], target['DEC'], target['Program'], target['Ebmv'],
                                       target['maxLen'], target['MoonFrac'], target['DESsn2'], target['Status'],
                                       target['Exposure'], target['obsSN2'], tbase))
                    else:
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
                        prihdr['EXPTIME '] = target['Exposure']
                        prihdr['OBSSN2  '] = target['obsSN2']
                        t = Time(mjd, format = 'mjd')
                        tbase = str(t.isot)
                        nt = len(tbase)
                        prihdr['DATE-OBS'] = tbase
                        filename = day_stats['dirName'] + '/desi-exp-' + ocnt.update() + '.fits'
                        prihdu = pyfits.PrimaryHDU(header=prihdr)
                        prihdu.writeto(filename, clobber=True)
                else:
                    # Try another target?
                    # Observe longer split into modulo(max_len)
                    mjd += 0.25/24.0
            else:
                mjd += 0.25/24.0
            # Check time
            if mjd > day_stats['MJDsunrise']:
                nightOver = True

    if tableOutput and len(obsList) > 0:
        filename = 'obslist' + day_stats['dirName'] + '.fits'
        cols = np.rec.array(obsList,
                           names = ('TILEID  ',
                                    'RA      ',
                                    'DEC     ',
                                    'PROGRAM ',
                                    'EBMV    ',
                                    'MAXLEN  ',
                                    'MOONFRAC',
                                    'DESSN2  ',
                                    'STATUS  ',
                                    'EXPTIME ',
                                    'OBSSN2  ',
                                    'DATE-OBS'),
                            formats = ['i4', 'f8', 'f8', 'a8', 'f8', 'f8', 'f8', 'f8', 'i4', 'f8', 'f8', 'a24'])
        tbhdu = pyfits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(filename, clobber=True)
    
    return tilesObserved

def surveySim(startday, startmonth, startyear, endday, endmonth, endyear):

    sp = surveyPlan()
    day0 = Time(datetime(startyear, startmonth, startday, 12, 0, 0))
    mjd_start = day0.mjd
    w = weatherModule(day0)
    ocnt = obsCount()

    tile_file = 'tiles_observed.fits'
    if os.path.exists(tile_file):
        tilesObserved = Table.read(tile_file, format='fits')
    else:
        print("The survey will start from scratch.")
        tilesObserved = Table(names=('TILEID', 'STATUS'), dtype=('i8', 'i4'))
        tilesObserved.meta['MJDBEGIN'] = mjd_start
    print ("Everything initialised.")

    cal = obsCalendar(startday, startmonth, startyear, endday, endmonth, endyear)
    for day in cal:
        ntodate = len(tilesObserved)
        t = Time(day['MJDsunset'], format = 'mjd')
        w.resetDome(t)
        tiles_todo, obsplan = sp.afternoonPlan(day, tilesObserved)
        tilesObserved = nightOps(day, obsplan, w, ocnt, tilesObserved)
        t = Time(day['MJDsunset'], format = 'mjd')
        print ('On the night starting ', t.iso, ', we observed ', len(tilesObserved)-ntodate, ' tiles.')
        if tiles_todo == 0:
            break

    tilesObserved.write(tile_file, format='fits', overwrite=True)
