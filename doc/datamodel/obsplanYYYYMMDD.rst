===============
obsplanYYYYMMDD
===============

:Summary: Contains the afternoon plan.
:Naming Convention: ``obsplanYYYYMMDD.fits``, where YYYYMMDD is the date of the
    start of the night of observations.
:Regex: ``obsplan20([0-9]{2})(0[1-9]|1[0-2])(0[1-9]|1[0-9]|2[0-9]|3[0-1])\.fits``
:File Type: FITS, 8-560 KB

Contents
========

====== ======= ======== ===================
Number EXTNAME Type     Contents
====== ======= ======== ===================
HDU0_          IMAGE    N/A
HDU1_          BINTABLE List of tiles in plan
====== ======= ======== ===================


FITS Header Units
=================

HDU0
----

N/A.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

======== ================= ===== =======
KEY      Example Value     Type  Comment
======== ================= ===== =======
MOONFRAC 0.144870860633956 float Moon illumination fraction
======== ================= ===== =======

Empty HDU.

HDU1
----

List of tiles to observe in order of priority.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

====== ============= ==== =====================
KEY    Example Value Type Comment
====== ============= ==== =====================
NAXIS1 56            int  length of dimension 1
NAXIS2 9607          int  length of dimension 2
====== ============= ==== =====================

Required Data Table Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

============= ======= ===== ===========
Name          Type    Units Description
============= ======= ===== ===========
TILEID        int32         DESI tile ID
RA            float64 deg   Right ascencion
DEC           float64 deg   Declination
EBV_MED       float64       Extinction
LSTMIN        float32 deg   Start of LST observing window
LSTMAX        float32 deg   End of LST observing window
MAXEXPLEN     float32 sec   Maximum exposure lenght
PRIORITY      int32         Between 0(high priority) - 10(low priority)
STATUS        int32         0 (unobserved) or 1 (partially observed)
PROGRAM       char[6]       'DARK', 'BRIGHT' or 'GRAY'
OBSCONDITIONS int16         0 (DARK), 1 (GRAY), 2 (BRIGHT) or 3 (POOR)
============= ======= ===== ===========


Notes and Examples
==================

Status==2 tiles (observed) should not be in this list.

