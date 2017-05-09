# surveysim

This package simulates the execution of the DESI main survey (ELG, LRG, QSO) and BGS survey.  For example, to run the first year of the nominal survey:

    % surveysim --start 2019-08-28 --end 2020-07-13 --verbose

For details on all the available options use:

    % surveysim --help

A plotting utility is provided to look at the progression of the survey and various metrics.  To run the plotting tool:

	>>> from surveysim.plotsurvey import plotsurvey
	>>> plotsurvey("obslist{_all|YYYYMMDD}.fits", plot_type='t', program='m')

The default filename is ./obslist_all.fits; plot_type is either 'f' (footprint, default), 'h' (histograms), 't' (time evolution) or 'e' (exposure time); and program us either 'm' (main survey), 'b' (BGS) or 'a' (all).
