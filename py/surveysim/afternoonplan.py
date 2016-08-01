import numpy as np
import astropy.io.fits as pyfits
from astropy.time import Time
from pkg_resources import resource_filename

Lat_KPNO_deg = 31.0 + (57.0 + 50.3/60.0)/60.0

class surveyPlan:

    def __init__(self):
    # Read list of DESI tiles
    # Format of prototype file:
    # name:
    #     ['TILEID', 'RA', 'DEC', 'PASS', 'IN_DESI', 'EBV_MED', 'AIRMASS', 'EXPOSEFAC', 'STAR_DENSITY', 'PROGRAM', 'OBSCONDITIONS']
    # format:
    #     ['J', 'D', 'D', 'I', 'I', 'E', 'E', 'E', 'E', '6A', 'I']
    # unit:
    #     ['', '', '', '', '', '', '', '', '']
    # null:
    #     ['', '', '', '', '', '', '', '', '']
    # bscale:
    #     ['', '', '', '', '', '', '', '', '']
    # bzero:
    #     ['', '', '', '', '', '', '', '', '']
    # disp:
    #     ['', '', '', '', '', '', '', '', '']
    # start:
    #     ['', '', '', '', '', '', '', '', '']
    # dim:
    #     ['', '', '', '', '', '', '', '', '']
        hdulist0 = pyfits.open(resource_filename('surveysim', 'data/desi-tiles.fits'))
        tiledata = hdulist0[1].data
        tileID = tiledata.field('TILEID')
        RA = tiledata.field('RA')
        DEC = tiledata.field('DEC')
        Pass = tiledata.field('PASS')
        InDESI = tiledata.field('IN_DESI')
        Ebmv = tiledata.field('EBV_MED')
        AM = tiledata.field('AIRMASS')
        expFac = tiledata.field('EXPOSEFAC')
        starDensity = tiledata.field('STAR_DENSITY')
        program = tiledata.field('PROGRAM')
        obsconds = tiledata.field('OBSCONDITIONS')
        obstime = tiledata.field('OBSTIME')
        HA = tiledata.field('HA')
        hdulist0.close()

        self.tileID = tileID.compress((InDESI==1).flat) #Assuming 0=out, 1=in
        self.RA = RA.compress((InDESI==1).flat)
        self.DEC = DEC.compress((InDESI==1).flat)
        self.Pass = Pass.compress((InDESI==1).flat)
        self.Ebmv = Ebmv.compress((InDESI==1).flat)
        self.maxExpLen = 2.0 * obstime.compress((InDESI==1).flat) # This assumes bright time program
        self.starDensity = starDensity.compress((InDESI==1).flat)
        self.program = program.compress((InDESI==1).flat)
        self.obsconds = obsconds.compress((InDESI==1).flat)
        LST = self.RA + HA.compress((InDESI==1).flat)
        self.LSTmin = LST - 15.0
        for i in range(len(self.LSTmin)):
            if self.LSTmin[i] < 0.0:
                self.LSTmin[i] += 360.0
            elif self.LSTmin[i] >360.0:
                self.LSTmin[i] -= 360.0
        self.LSTmax = LST + 15.0
        for i in range(len(self.LSTmax)):
            if self.LSTmax[i] < 0.0:
                self.LSTmax[i] += 360.0
            elif self.LSTmax[i] > 360.0:
                self.LSTmax[i] -= 360.0

        self.status = np.zeros(len(self.tileID))
        self.priority = np.ones(len(self.tileID))

    def afternoonPlan(self, day_stats, tiles_observed):
        """
        All the file names are hard coded, so there is no need to
        have them as arguments to this function.  The output file
        name has format obsplanYYYYMMDD.fits .
        """

        nto = len(tiles_observed)
        for i in range(nto):
            j = np.where(self.tileID == tiles_observed['TILEID'][i])
            self.status[j] = tiles_observed['STATUS'][i]
        
        # table: ra, dec, tileid, priority, design HA, estimated exp time, LST begin/end and or min/max, templates, programme name, S/N2
        # FOR NOW RETURN ALL UNOBSERVED TILES
        a1 = self.tileID.compress((self.status<2).flat)
        col1 = pyfits.Column(name='TILEID', format='J', array=a1)

        a2 = self.RA.compress((self.status<2).flat)
        col2 = pyfits.Column(name='RA', format='D', array=a2)

        a3 = self.DEC.compress((self.status<2).flat)
        col3 = pyfits.Column(name='DEC', format='D', array=a3)

        a4 = self.Ebmv.compress((self.status<2).flat)
        col4 = pyfits.Column(name='EBV_MED', format='E', array=a4)

        a5 = self.LSTmin.compress((self.status<2).flat)
        col5 = pyfits.Column(name='LSTMIN', format='E', array=a5)

        a6 = self.LSTmax.compress((self.status<2).flat)
        col6 = pyfits.Column(name='LSTMAX', format='E', array=a6)

        a7 = self.maxExpLen.compress((self.status<2).flat)
        col7 = pyfits.Column(name='MAXEXPLEN', format='E', array=a7)

        a8 = self.priority.compress((self.status<2).flat)
        col8 = pyfits.Column(name='PRIORITY', format='J', array=a8)
        
        a9 = self.status.compress((self.status<2).flat)
        col9 = pyfits.Column(name='STATUS', format='J', array=a9)

        a10 = self.program.compress((self.status<2).flat)
        col10 = pyfits.Column(name='PROGRAM', format='6A', array=a10)

        a11 = self.obsconds.compress((self.status<2).flat)
        col11 = pyfits.Column(name='OBSCONDITIONS', format='I', array=a11)

        cols = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11])
        tbhdu = pyfits.BinTableHDU.from_columns(cols)

        prihdr = pyfits.Header()
        prihdr['MOONFRAC'] = day_stats['MoonFrac']
        prihdr['MOONRA  '] = day_stats['MoonRA']
        prihdr['MOONDEC '] = day_stats['MoonDEC']
        prihdu = pyfits.PrimaryHDU(header=prihdr)

        t = Time(day_stats['MJDsunset'], format = 'mjd')
        utc = t.iso
        filename = 'obsplan' + day_stats['dirName'] + '.fits'
        thdulist = pyfits.HDUList([prihdu, tbhdu])
        thdulist.writeto(filename, clobber=True)

        tilesTODO = len(a9)

        return tilesTODO, filename

        
    
