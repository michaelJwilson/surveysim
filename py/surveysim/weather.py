from __future__ import print_function, division, absolute_import

from datetime import datetime

import numpy as np

import astropy.time
import astropy.table
import astropy.units as u

import desiutil.log

import desimodel.seeing

import desisurvey.config
import desisurvey.utils


# Percentage probability of the dome being closed due to weather
# during each calendar month.
dome_closed_probability = [
    35.24, 44.14, 27.68, 26.73, 14.22, 15.78,
    55.92, 48.75, 29.45, 24.44, 24.86, 34.74]


def sample_gaussian_random_process(
    n_sample, dt_sec=300., mean= 0., sigma=1., tau1=0.025, tau2=5., gen=None):
    """Generate samples of a Gaussian random process.

    Samples are correlated in time according to the value of two
    time constants.

    This is a simplified version of desimodel.seeing.sample() and should
    probably go into desiutils somewhere.

    Parameters
    ----------
    n_sample : int
        Number of samples to generate. Must be at least 2.
    dt_sec : float
        Elapsed time in seconds between generated samples.
    mean : float
        Mean value of the generated samples.
    sigma : float
        Standard deviation of the generated samples.
    tau1 : float
        Time constant associated with the autocorrelation power spectral
        density, in units of days.
    tau2 : float
        Time constant associated with the autocorrelation power spectral
        density, in units of days.
    gen : numpy.random.RandomState
        Random number generator to use for reproducible samples.

    Returns
    -------
    array:
        Array of sampled values on a uniform time grid.
    """
    if n_sample < 2:
        raise ValueError('n_sample must be at least 2.')

    # Build a linear grid of frequencies present in the Fourier transform
    # of the requested time series.  Frequency units are 1/day.
    dt_day = dt_sec / (24. * 3600.)
    df_day = 1. / (n_sample * dt_day)
    f_grid = np.arange(1 + (n_sample // 2)) * df_day

    # Tabulate the power spectral density at each frequency.
    omega = 2 * np.pi * f_grid
    psd = 1. / (1. + omega * tau1) / (1. + omega * tau2)
    # Force the mean to zero.
    psd[0] = 0.
    # Force the variance to one.
    psd[1:] /= psd[1:].sum() * df_day ** 2

    # Sample the Gaussian random process.
    n_psd = len(psd)
    x_fft = np.ones(n_psd, dtype=complex)
    x_fft[1:-1].real = gen.normal(size=n_psd - 2)
    x_fft[1:-1].imag = gen.normal(size=n_psd - 2)
    x_fft *= np.sqrt(psd) / (2 * dt_day)
    x_fft[0] *= np.sqrt(2)
    x = np.fft.irfft(x_fft, n_sample)

    # Transform samples to the requested mean and sigma.
    x = mean + sigma * x

    return x


class Weather(object):
    """Simulate weather conditions affecting observations.

    Seeing and transparency values are stored with 32-bit floats to save
    some memory.

    Parameters
    ----------
    start_date : datetime.date or None
        Survey starts on the evening of this date. Use the ``first_day``
        config parameter if None (the default).
    stop_date : datetime.date or None
        Survey stops on the morning of this date. Use the ``last_day``
        config parameter if None (the default).
    time_step : astropy.units.Quantity
        Time step for calculating updates. Must evenly divide 24 hours.
    gen : numpy.random.RandomState or None
        Random number generator to use for reproducible samples. Will be
        initialized (un-reproducibly) if None.
    """
    def __init__(self, start_date=None, stop_date=None, time_step=5 * u.min,
                 gen=None):
        self.log = desiutil.log.get_logger()
        config = desisurvey.config.Configuration()

        if gen is None:
            self.log.warn('Will generate unreproducible random numbers.')
            gen = np.random.RandomState()

        # Use our config to set any unspecified dates.
        if start_date is None:
            start_date = config.first_day()
        if stop_date is None:
            stop_date = config.last_day()
        num_nights = (stop_date - start_date).days
        if num_nights <= 0:
            raise ValueError('Expected start_date < stop_date.')

        # Check that the time step evenly divides 24 hours.
        steps_per_day = int(round((1 * u.day / time_step).to(1).value))
        if not np.allclose((steps_per_day * time_step).to(u.day).value, 1.):
            raise ValueError(
                'Requested time_step does not evenly divide 24 hours: {0}.'
                .format(time_step))

        # Calculate the number of times where we will tabulate the weather.
        num_rows = num_nights * steps_per_day
        self._table = astropy.table.Table()

        # Initialize column of MJD timestamps.
        t0 = desisurvey.utils.local_noon_on_date(start_date)
        times = t0 + (np.arange(num_rows) / float(steps_per_day)) * u.day
        self._table['mjd'] = times.mjd

        # Decide whether the dome is opened on each night.
        # We currently assume this is fixed for a whole night, but
        # tabulate the status at each time so that this could be
        # updated in future to simulate partial-night weather outages.
        self._table['open'] = np.ones(num_rows, bool)
        for i in range(num_nights):
            ij = i * steps_per_day
            month = times[ij].datetime.month
            if 100 * gen.uniform() < dome_closed_probability[month - 1]:
                self._table['open'][ij:ij + steps_per_day] = False

        # Sample random seeing values.
        dt_days = 24 * 3600. / steps_per_day
        # It would be better if we could pass our gen to seeing.sample().
        seeing_seed = gen.randint(2 ** 31)
        self._table['seeing'] = desimodel.seeing.sample(
            num_rows, dt_days, seed=seeing_seed).astype(np.float32)

        # Sample transparency as the lognormal transform of a Gaussian
        # random process.  Mean and sigma are copied from the original code.
        self._table['transparency'] = np.clip(
            np.exp(sample_gaussian_random_process(
                num_rows, dt_days, mean=0.11111, sigma=0.33333, gen=gen)),
            0., 1.).astype(np.float32)

        self.start_date = start_date
        self.stop_date = stop_date
        self.num_nights = num_nights
        self.steps_per_day = steps_per_day

    def get(self, time):
        """Get the weather conditions at the specified time(s).

        Returns the conditions at the closest tabulated time, rather than
        using interpolation.

        Parameters
        ----------
        time : astropy.time.Time
            Time(s) when the simulated weather is requested.

        Returns
        -------
        table slice
            Slice of precomputed table containing row(s) corresponding
            to the requested time(s).
        """
        offset = np.floor(
            (time.mjd - self._table['mjd'][0]) * self.steps_per_day + 0.5
            ).astype(int)
        if np.any(offset < 0) or np.any(offset > len(self._table)):
            raise ValueError('Cannot get weather beyond tabulated range.')
        return self._table[offset]
