import numpy as np
import astropy.io.fits as pyfits
from astropy.time import Time

Lat_KPNO_deg = 31.0 + (57.0 + 50.3/60.0)/60.0

class surveyPlan:

    def __init__(self):
    # Read list of DESI tiles
    # Format of prototype file:
    # name:
    #     ['TILEID', 'RA', 'DEC', 'PASS', 'IN_DESI', 'EBV_MED', 'AIRMASS', 'EXPOSEFAC', 'STAR_DENSITY']
    # format:
    #     ['J', 'D', 'D', 'I', 'I', 'E', 'E', 'E', 'E']
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
        hdulist0 = pyfits.open('data/desi-tiles.fits')
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
        hdulist0.close()

        self.tileID = tileID.compress((InDESI==1).flat) #Assuming 0=out, 1=in
        self.RA = RA.compress((InDESI==1).flat)
        self.DEC = DEC.compress((InDESI==1).flat)
        self.Pass = Pass.compress((InDESI==1).flat)
        self.Ebmv = Ebmv.compress((InDESI==1).flat)
        self.maxExpLen = 1000.0 * expFac.compress((InDESI==1).flat) # This assumes bright time program
        self.starDensity = starDensity.compress((InDESI==1).flat)

        # Next versions of the file should have this instead of airmass
        phi = np.pi*Lat_KPNO_deg/180.0
        AM0 = AM.compress((InDESI==1).flat)
        cosHA = (1.0/AM0 - np.sin(phi)*np.sin(self.DEC)) / (np.cos(phi)*np.cos(self.DEC))
        iplus = np.where(cosHA > 1.0)
        cosHA[iplus] = 1.0
        iminus = np.where(cosHA < -1.0)
        cosHA[iminus] = -1.0
        LST = self.RA + np.arccos(cosHA)*180.0/np.pi
        self.LSTmin = LST - 15.0
        self.LSTmax = LST + 15.0

        self.status = np.empty(len(self.tileID))
        self.status.fill(-1)

        self.priority = np.ones(len(self.tileID))

    def afternoonPlan(self, day_stats, tiles_observed=None):
        """
        All the file names are hard coded, so there is no need to
        have them as arguments to this function.  The output file
        name has format obsplanYYYYMMDD.fits .
        """

        if tiles_observed != None:
            for i in xrange(len(tiles_observed)):
                j = np.where(self.tileID == tiles_observed[i]['tileID'])
                self.status[j] = 1
        
        # table: ra, dec, tileid, priority, design HA, estimated exp time, LST begin/end and or min/max, templates, programme name, S/N2
        # FOR NOW RETURN ALL UNOBSERVED TILES
        a1 = self.tileID.compress((self.status<0).flat)
        col1 = pyfits.Column(name='TILEID', format='J', array=a1)

        a2 = self.RA.compress((self.status<0).flat)
        col2 = pyfits.Column(name='RA', format='D', array=a2)

        a3 = self.DEC.compress((self.status<0).flat)
        col3 = pyfits.Column(name='DEC', format='D', array=a3)

        a4 = self.Ebmv.compress((self.status<0).flat)
        col4 = pyfits.Column(name='EBV_MED', format='E', array=a4)

        a5 = self.LSTmin.compress((self.status<0).flat)
        col5 = pyfits.Column(name='LSTMIN', format='E', array=a5)

        a6 = self.LSTmax.compress((self.status<0).flat)
        col6 = pyfits.Column(name='LSTMAX', format='E', array=a6)

        a7 = self.maxExpLen.compress((self.status<0).flat)
        col7 = pyfits.Column(name='MAXEXPLEN', format='E', array=a7)

        a8 = self.priority.compress((self.status<0).flat)
        col8 = pyfits.Column(name='PRIORITY', format='J', array=a8)
        
        a9 = self.status.compress((self.status<0).flat)
        col9 = pyfits.Column(name='STATUS', format='J', array=a9)

        cols = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])
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

        
    
