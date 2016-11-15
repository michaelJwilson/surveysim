==============
tiles_observed
==============

:Summary: Lists the tile ID of every tile observed so far.
:Naming Convention: ``tiles_observed.fits``
:Regex: ``tiles\_observed\.fits``
:File Type: FITS, 8 KB

Contents
========

====== ======= ======== ===================
Number EXTNAME Type     Contents
====== ======= ======== ===================
HDU0_          IMAGE    N/A
HDU1_          BINTABLE List of observed tiles
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

List of observed tiles.

Required Header Keywords
~~~~~~~~~~~~~~~~~~~~~~~~

======== ================= ===== =====================
KEY      Example Value     Type  Comment
======== ================= ===== =====================
NAXIS1   12                int   length of dimension 1
NAXIS2   141               int   length of dimension 2
MJDBEGIN 57749.79166666666 float Start day of the survey
======== ================= ===== =====================

Required Data Table Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

====== ===== ===== ===========
Name   Type  Units Description
====== ===== ===== ===========
TILEID int64       DESI tile ID
STATUS int32       1 (partially observed) or 2 (completed)
====== ===== ===== ===========


Notes and Examples
==================

