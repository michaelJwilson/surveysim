import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as ticker

"""
This function plots various quantities output from surveySim.
 -- Martin Landriau, August 2016
"""
def plotsurvey(filename='obslist_all.fits', plot_type='f'):
    t = Table.read(filename, format='fits')

    if plot_type == 'f':
        plt.figure(1)
        plt.subplot(111)
        ra = t['RA']
        dec = t['DEC']
        mjd = t['MJD']
        mjd_start = np.min(mjd)
        mjd -= mjd_start

        i1 = np.where( mjd/365.0 < 1.0 )
        i2 = np.where( (mjd/365.0 >= 1.0) & (mjd/365.0 < 2.0) )
        i3 = np.where( (mjd/365.0 >= 2.0) & (mjd/365.0 < 3.0) )
        i4 = np.where( (mjd/365.0 >= 3.0) & (mjd/365.0 < 4.0) )
        i5 = np.where( (mjd/365.0 >= 4.0) & (mjd/365.0 < 5.0) )
        y1 = plt.scatter(ra[i1], dec[i1], c='r')
        y2 = plt.scatter(ra[i2], dec[i2], c='b')
        y3 = plt.scatter(ra[i3], dec[i3], c='g')
        y4 = plt.scatter(ra[i4], dec[i4], c='y')
        y5 = plt.scatter(ra[i5], dec[i5], c='m')

        plt.xlabel('RA (deg)')
        plt.ylabel('DEC (deg)')
        plt.legend((y1, y2, y3, y4, y5), ('Year 1', 'Year 2', 'Year 3', 'Year 4', 'Year 5'), scatterpoints=1)

    elif plot_type == 'h':
        plt.figure(1)
        plt.subplot(231)
        x = t['EXPTIME']
        n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
        plt.xlabel('Exposure time (seconds)')
        plt.ylabel('Count')

        plt.subplot(232)
        x = t['SEEING']
        n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
        plt.xlabel('Seeing (arcseconds)')
        plt.ylabel('Count')

        plt.subplot(233)
        x = t['LINTRANS']
        n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
        plt.xlabel('Linear transparency')
        plt.ylabel('Count')

        plt.subplot(234)
        x = t['AIRMASS']
        n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
        plt.xlabel('Airmass')
        plt.ylabel('Count')

        plt.subplot(235)
        y = t['MOONALT']
        x = t['MOONFRAC'].compress((y>0.0).flat)
        n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
        plt.xlabel('Moon illumination fraction')
        plt.ylabel('Count')

        plt.subplot(236)
        y = t['MOONALT']
        x = t['MOONDIST'].compress((y>0.0).flat)
        n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
        plt.xlabel('Distance from the Moon (deg)')
        plt.ylabel('Count')

    elif plot_type == 't':
        mjd = t['MJD']
        plt.figure(1)

        plt.subplot(221)
        y = t['MOONALT']
        plt.plot(mjd, y, linestyle='-', color='black')
        plt.xlabel('MJD')
        plt.ylabel('Moon elevation (degrees)')

        plt.subplot(222)
        y = t['SEEING']
        plt.plot(mjd, y, linestyle='-', color='black')
        plt.xlabel('MJD')
        plt.ylabel('Seeing (arcseconds)')

        plt.subplot(223)
        y = t['LINTRANS']
        plt.plot(mjd, y, linestyle='-', color='black')
        plt.xlabel('MJD')
        plt.ylabel('Linear transparency')

        plt.subplot(224)
        y = np.arange(len(mjd)) + 1
        plt.plot(mjd, y, linestyle='-', color='black')
        plt.xlabel('MJD')
        plt.ylabel('Number of tiles observed')



    plt.show()

