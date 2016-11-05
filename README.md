# surveysim

The full path to surveysim/py should be added to PYTHONPATH.

To run, for example:

	>>> from surveysim.wrapper import surveySim
	>>> surveySim((2016, 12, 27), (2017, 1, 4), 123456)

The optional seed for the weather module's random number generator has to
be an int or array_like convertible to an unsigned 32-bit integer.

To run the plotting tool:

	>>> from surveysim.plotsurvey import plotsurvey
	>>> plotsurvey("your_output_file.fits, plot_type='t')

The default filename is ./obslist_all.fits and plot_type is either 'f' (footprint, default), 'h' (histograms) or 't' (time evolution).

