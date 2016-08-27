import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as ticker

"""
This is a short script to plot the evolution of the survey
from the file tiles_observed.dat (same as tiles_observed.fits,
but in ASCII format) output from surveySim.
 -- Martin Landriau, August 2016
"""

t = Table.read('obslist_all.fits', format='fits')

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
x = t['MOONFRAC']
n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
plt.xlabel('Moon illumination fraction')
plt.ylabel('Count')

plt.subplot(326)
x = t['MOONDIST']
n, bins, patches = plt.hist(x, 20, facecolor='0.5', alpha=0.75)
plt.xlabel('Distance from the Moon (deg)')
plt.ylabel('Count')

plt.show()

