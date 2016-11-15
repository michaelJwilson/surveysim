===========
obslist_all
===========

:Summary: List of tiles observed, along with tile info and observing conditions.
:Naming Convention: ``obslist_all.fits``
:Regex: ``obslist\_all\.fits``
:File Type: FITS, 1 MB 

Contents
========

====== ======= ======== ===================
Number EXTNAME Type     Contents
====== ======= ======== ===================
HDU0_          IMAGE    N/A
HDU1_          BINTABLE Observed tiles
====== ======= ======== ===================


FITS Header Units
=================

HDU0
----

N/A.

This HDU has no non-standard required keywords.

Empty HDU.

HDU1
----

Table of observed tiles.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

====== ============= ==== =====================
KEY    Example Value Type Comment
====== ============= ==== =====================
NAXIS1 152           int  length of dimension 1
NAXIS2 10092         int  length of dimension 2
====== ============= ==== =====================

Required Data Table Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

======== ======== ===== ===========
Name     Type     Units Description
======== ======== ===== ===========
TILEID   int32          DESI tile ID
RA       float64  deg   right ascension
DEC      float64  deg   declination
PROGRAM  char[8]        'DARK', 'BRIGHT' or 'GRAY'
EBMV     float64        Extinction
MAXLEN   float64  sec   Maximum exposure length
MOONFRAC float64        Fraction of the Moon's illuminated surface
MOONDIST float64  deg   Distance between the Moon and the tile centre.
MOONALT  float64  deg   Elevation of the Moon
SEEING   float64  "     Seeing
LINTRANS float64        Linear transparency
AIRMASS  float64        Airmass
DESSN2   float64        Design (S/N)^2
STATUS   int32          2 (fully observed)
EXPTIME  float64  sec   Exposure time
OBSSN2   float64        Achieved (S/N)^2
DATE-OBS char[24]       Date and time of observation
MJD      float64  days  Modified Julian Date of observation 
======== ======== ===== ===========


Notes and Examples
==================

Tiles with status==0 (unobserved) or status==1 (partially observed) should not
be in this list.

