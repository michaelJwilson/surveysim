import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as ticker

"""
This function plots various quantities output from surveySim.
 -- Martin Landriau, August 2016
"""
def plotsurvey(filename='obslist_all.fits'):
    t = Table.read(filename, format='fits')

    plt.figure(1)
    plt.subplot(321)
    plt.scatter(t['RA'], t['DEC'])
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')

    plt.subplot(322)
    x = t['SEEING']
    n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
    plt.xlabel('Seeing (arcseconds)')
    plt.ylabel('Count')

    plt.subplot(323)
    x = t['LINTRANS']
    n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
    plt.xlabel('Linear transparency')
    plt.ylabel('Count')

    plt.subplot(324)
    x = t['AIRMASS']
    n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
    plt.xlabel('Airmass')
    plt.ylabel('Count')

    plt.subplot(325)
    y = t['MOONALT']
    x = t['MOONFRAC'].compress((y>0.0).flat)
    n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
    plt.xlabel('Moon illumination fraction')
    plt.ylabel('Count')

    plt.subplot(326)
    y = t['MOONALT']
    x = t['MOONDIST'].compress((y>0.0).flat)
    n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
    plt.xlabel('Distance from the Moon (deg)')
    plt.ylabel('Count')

    plt.show()

