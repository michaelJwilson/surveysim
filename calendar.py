from astropy.time import Time, TimezoneInfo
from datetime import datetime
import astropy.units as u
import numpy as np

"""
Line format
0, 1, 2: evening day, month, date
3: /
4, 5, 6: morning day, month, date
7: MJD of local midnight
8, 9, 10: LMST of midnight hours, minutes, seconds
11, 12: Sun set hours, minutes
13, 14: evening twilight end hours, minutes
15, 16: morning twilight begin hours, minutes
17, 18: Sun rise hours, minutes
19, 20: evening twilight LST hours, minutes
21, 22: morning twilight LST hours, minutes
23, 24: Moon rise hours, minutes                   23: ....                                  23, 24: Moon rise hours, minutes
25, 26: Moon set hours, minutes                    24, 25: Moon set hours, minutes           25: ....
27: Moon illumination fraction                     26: Moon illumination fraction            26: Moon illumination fraction
28, 29: Moon RA hours, minutes                     27, 28: Moon RA hours, minutes            27, 28: Moon RA hours, minutes
30: -
32, 33: Moon DEC degrees, minutes                  31, 32: Moon DEC degrees, minutes         31, 32: Moon DEC degrees, minutes
"""
daymonth = np.array([31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) # This is only valid for leap years like 2016.

def monthNum(monthname):
    month = 0
    if monthname == 'Jan': month = 1
    elif monthname == 'Feb': month = 2
    elif monthname == 'Mar': month = 3
    elif monthname == 'Apr': month = 4
    elif monthname == 'May': month = 5
    elif monthname == 'Jun': month = 6
    elif monthname == 'Jul': month = 7
    elif monthname == 'Aug': month = 8
    elif monthname == 'Sep': month = 9
    elif monthname == 'Oct': month = 10
    elif monthname == 'Nov': month = 11
    elif monthname == 'Dec': month = 12
    return month

def obsCalendar(calfile):
    utc_minus7 = TimezoneInfo(utc_offset=-7*u.hour)
    year = 2016 # Testing calendar with only one year
    cal = []
    with open(calfile,'r') as f:
        while True:
            line = f.readline()
            if not line: break
            if line[0] == '#': continue
            entries = line.split()
            month = monthNum(entries[1])
            day = int(entries[2])
            MJDmidnight = float(entries[7])
            sunset_h = int(entries[11])
            sunset_m = int(entries[12])
            Tsunset = Time(datetime(year, month, day, sunset_h, sunset_m, 0, tzinfo=utc_minus7))
            etwi_h = int(entries[13])
            etwi_m = int(entries[14])
            Tetwi = Time(datetime(year, month, day, etwi_h, etwi_m, 0, tzinfo=utc_minus7))
            mtwi_h = int(entries[15])
            mtwi_m = int(entries[16])
            yearnew = year
            monthnew = month
            daynew = day + 1
            if daynew > daymonth[month-1]:
                daynew = 1
                monthnew = month + 1
                if monthnew > 12:
                    monthnew = 1
                    yearnew = year +1
            Tmtwi = Time(datetime(yearnew, monthnew, daynew, mtwi_h, mtwi_m, 0, tzinfo=utc_minus7))
            sunrise_h = int(entries[17])
            sunrise_m = int(entries[18])
            Tsunrise = Time(datetime(yearnew, monthnew, daynew, sunrise_h, sunrise_m, 0, tzinfo=utc_minus7))
            if entries[23] == ".....":
                MJDmoonrise = -1.0
                moonset_h = int(entries[24])
                moonset_m = int(entries[25])
                yearnew = year
                monthnew = month
                daynew = day
                if moonset_h < 12:
                    daynew = day + 1
                    if daynew > daymonth[month-1]:
                        daynew = 1
                        monthnew = month + 1
                        if monthnew > 12:
                            monthnew = 1
                            yearnew = year +1
                MJDmoonset = (Time(datetime(yearnew, monthnew, daynew, moonset_h, moonset_m, 0, tzinfo=utc_minus7))).mjd
            else:
                moonrise_h = int(entries[23])
                moonrise_m = int(entries[24])
                yearnew = year
                monthnew = month
                daynew = day
                if moonrise_h < 12:
                    daynew = day + 1
                    if daynew > daymonth[month-1]:
                        daynew = 1
                        monthnew = month + 1
                        if monthnew > 12:
                            monthnew = 1
                            yearnew = year +1
                MJDmoonrise = (Time(datetime(yearnew, monthnew, daynew, moonrise_h, moonrise_m, 0, tzinfo=utc_minus7))).mjd
                if entries[25] == ".....":
                    MJDmoonset = -1.0
                else:
                    moonset_h = int(entries[25])
                    moonset_m = int(entries[26])
                    yearnew = year
                    monthnew = month
                    daynew = day
                    if moonset_h < 12:
                        daynew = day + 1
                        if daynew > daymonth[month-1]:
                            daynew = 1
                            monthnew = month + 1
                            if monthnew > 12:
                                monthnew = 1
                                yearnew = year +1
                    MJDmoonset = (Time(datetime(yearnew, monthnew, daynew, moonset_h, moonset_m, 0, tzinfo=utc_minus7))).mjd
            if entries[23] == "....." or entries[25] == ".....":
                moonIllum = float(entries[26])
                moonRA = float(entries[27]) + float(entries[28])/60.0
                moonRA = moonRA * 15.0
                if entries[29] == "-":
                    moonDEC = -1.0 * (float(entries[30]) + float(entries[31])/100.0)
                else:
                    moonDec = float(entries[29]) + float(entries[30])/100.0
            else:
                moonIllum = float(entries[27])
                moonRA = float(entries[28]) + float(entries[29])/60.0
                if entries[30] == "-":
                    moonDEC = -1.0 * (float(entries[31]) + float(entries[32])/100.0)
                else:
                    moonDec = float(entries[30]) + float(entries[31])/100.0
            cal.append( {'MJDmidnight': MJDmidnight,
                         'MJDsunset': Tsunset.mjd,
                         'MJDsunrise': Tsunrise.mjd,
                         'MJDetwi': Tetwi.mjd,
                         'MJDmtwi': Tmtwi.mjd,
                         'MJDmoonrise': MJDmoonrise,
                         'MJDmoonset': MJDmoonset,
                         'MoonFrac': moonIllum,
                         'MoonRA': moonRA,
                         'MoonDEC': moonDEC} )
    return cal
