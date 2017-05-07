====================
surveysim change log
====================

0.5.0 (unreleased)
------------------

* Use desisurvey.config to manage all non-simulation configuration data.
* Unify different output files with overlapping contents into single output
  managed by desisurvey.progress.
* Overhaul of weather simulator to generate continuous stationary time series
  that are independent of the observing sequence.  Use desimodel.seeing.
* Simulate multiple exposures for cosmics and more realistic overhead.
* Clean up of README, docstrings, imports, unit tests, requirements, unused code.

0.4.1 (2017-04-13)
------------------

* Fixed package names to work with desisurvey >= 0.4.0

0.4.0 (2017-04-04)
------------------

* Adds unit tests
* removes data/tile-info.fits (not used here; was moved to desisurvey)
* adds nightops.py (from desisurvey, used here but not there)
* create surveysim command-line script
* use new desisurvey config machinery (first steps, in progress)

0.3.1 (2016-12-21)
------------------

* Fixed outlier HA tile assignments around RA 200-220 (PR #26)
* Added 7 day shutdown around full moon (PR #25)

0.3.0 (2016-11-29)
------------------

* Moved non-simulation specific parts to desisurvey

0.2.0 (2016-11-18)
------------------

* Modified some file names
* Moved some functions from one file to another

0.1.1 (2016-11-14)
------------------

* fixed crash at end and data/ install (PR #3)
* initial tests for NERSC install

0.1.0 and prior
---------------

* No changes.rst yet
