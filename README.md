# surveysim

The full path to surveysim/py should be added to PYTHONPATH.

To run, for example:

	>>> from surveysim.surveysim import surveySim
	>>> surveySim((2016, 12, 27), (2017, 1, 4), seed=123456, use_jpl=False)

The optional seed for the weather module's random number generator has to
be an int or array_like convertible to an unsigned 32-bit integer; use_jpl specifies which version of avoidobject.py to use.

