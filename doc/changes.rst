====================
surveysim change log
====================

0.10.0 (unreleased)
-------------------

* Refactor desisurvey.schedule -> desisurvey.old.schedule
* Add new modules: exposures, stats.
* Requires desisurvey 0.11.0.

0.9.2 (2018-10-02)
------------------

* Replay historical Mayall daily weather.
* Implement partial-night dome closure.
* Requires desimodel >= 0.9.8 and desisurvey >= 0.10.4.

0.9.1 (2018-06-27)
------------------

* Do arc exposures before flat exposures (PR `#57`_).

.. _`#57`: https://github.com/desihub/surveysim/pull/57

0.9.0 (2017-11-09)
------------------

* Add ``surveysim.util.add_calibration_exposures()``, to add simulated
  calibration exposures to a set of science exposures (PR `#55`_).

.. _`#55`: https://github.com/desihub/surveysim/pull/55

0.8.2 (2017-10-09)
------------------

* Use new desisurvey config api (requires desisurvey >= 0.9.3)
* Add support for optional depth-first survey strategy.
* Docs now auto-generated at http://surveysim.readthedocs.io/en/latest/

0.8.1 (2017-09-20)
------------------

* Adds surveysim --config-file option (PR `#49`_); requires desisurvey/0.9.1.

.. _`#49`: https://github.com/desihub/surveysim/pull/49

0.8.0 (2017-09-11)
------------------

* Track API changes in desisurvey 0.9.0.
* The surveysim script is now called once per night, alternating with a
  surveyplan script that lives in the desisurvey package.
* See https://www.youtube.com/watch?v=vO1QZD_aCIo for a visualization of the
  full 5-year survey simulation that matches DESI-doc-1767-v3.

0.7.1 (2017-08-07)
------------------

* Use new desimodel.weather to randomly sample seeing and transparency.
  Requires desimodel >= 0.8.0.

0.7.0 (2017-06-18)
------------------

* First implementation of fiber-assignment groups and priorities.
* Integration with the new desisurvey surveyplan script.
* Create tutorial document and sample automation scripts.

0.6.0 (2017-06-05)
------------------

* Add strategy, weights options to surveysim script.
* Add hooks for using greedy scheduler
* Terminate exposures at sunset

0.5.0 (2017-05-10)
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
