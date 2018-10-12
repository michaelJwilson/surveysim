"""Record simulated nightly statistics by pass.
"""
from __future__ import print_function, division, absolute_import

import numpy as np

import desisurvey.config
import desisurvey.utils

class SurveyStatistics(object):
    """Collect nightly statistics by pass.
    """
    def __init__(self, scheduler):
        self.scheduler = scheduler
        config = desisurvey.config.Configuration()
        self.start_date = config.first_day()
        self.stop_date = config.last_day()
        self.num_nights = (self.stop_date - self.start_date).days
        # Build our internal array.
        dtype = []
        for name in 'MJD', 'tsched',:
            dtype.append((name, np.float))
        nprograms = len(scheduler.program_passes)
        for name in 'topen', 'tdead',:
            dtype.append((name, np.float, (nprograms,)))
        npasses = len(scheduler.pass_program)
        for name in 'tscience', 'tsetup', 'tsplit',:
            dtype.append((name, np.float, (npasses,)))
        for name in 'completed', 'nexp', 'nsetup', 'nsplit', 'nsetup_abort', 'nsplit_abort',:
            dtype.append((name, np.int32, (npasses,)))
        self._data = np.empty(self.num_nights, dtype)
        self.reset()

    def reset(self):
        self._data[:] = 0
        # Initialize local-noon MJD timestamp for each night.
        first_noon = desisurvey.utils.local_noon_on_date(self.start_date).mjd
        self._data['MJD'] = first_noon + np.arange(self.num_nights)
        
    @property
    def nexp(self):
        return self._data['nexp'].sum()
        
    def get_night(self, night):
        night = desisurvey.utils.get_date(night)
        assert night < self.stop_date
        idx = (night - self.start_date).days
        return self._data[idx]
    
    def validate(self):
        D = self._data
        # Every exposure must be preceded by a setup or split.
        if not np.all(D['nexp'] == D['nsplit'] + D['nsetup']):
            return False
        # Sum live time per program over nights.
        tlive = (D['topen'] - D['tdead']).sum(axis=1)
        # Sum time spent in each state per pass over nights.
        ttotal = (D['tsetup'] + D['tscience'] + D['tsplit']).sum(axis=1)
        return np.allclose(tlive, ttotal)

    def summarize(self):
        assert self.validate()
        D = self._data
        tsched = 24 * D['tsched'].sum()
        topen = 24 * D['topen'].sum()
        tscience = 24 * D['tscience'].sum()
        print('Scheduled {:.3f} hr Open {:.3f}% Live {:.3f}%'.format(
            tsched, 100 * topen / tsched, 100 * tscience / topen))
        print('=' * 82)
        print('PROG PASS    TILES  NEXP SETUP ABT SPLIT ABT    TEXP TSETUP TSPLIT   TOPEN  TDEAD')
        print('=' * 82)
        # Summarize by pass.
        for pidx, program in enumerate(self.scheduler.programs):
            ntiles_p, ndone_p, nexp_p, nsetup_p, nsplit_p, nsetup_abort_p, nsplit_abort_p = [0] * 7
            tscience_p, tsetup_p, tsplit_p = [0.] * 3
            passes = []
            ntiles_all = 0
            for passnum in list(self.scheduler.program_passes[program]) + [' ']:
                if passnum == ' ':
                    sel = passes
                    ntiles = ntiles_all
                else:
                    sel = passnum
                    ntiles = np.count_nonzero(self.scheduler.passnum == passnum)
                    ntiles_all += ntiles
                    passes.append(passnum)
                ndone = D['completed'][:, sel].sum()
                nexp = D['nexp'][:, sel].sum()
                nsetup = D['nsetup'][:, sel].sum()
                nsplit = D['nsplit'][:, sel].sum()
                nsetup_abort = D['nsetup_abort'][:, sel].sum()
                nsplit_abort = D['nsplit_abort'][:, sel].sum()
                tscience = 86400 * D['tscience'][:, sel].sum() / ndone
                tsetup = 86400 * D['tsetup'][:, sel].sum() / ndone
                tsplit = 86400 * D['tsplit'][:, sel].sum() / ndone
                line = '{:6s} {} {:4d}/{:4d} {:5d} {:5d} {:3d} {:5d} {:3d} {:6.1f}s {:5.1f}s {:5.1f}s'.format(
                    program, passnum, ndone, ntiles, nexp, nsetup, nsetup_abort, nsplit, nsplit_abort, tscience, tsetup, tsplit)
                if passnum == ' ':
                    # Open and deadtime are accumulated by program, not pass.
                    topen = 86400 * D['topen'][:, pidx].sum() / ndone
                    tdead = 86400 * D['tdead'][:, pidx].sum() / ndone
                    line += ' {:6.1f}s {:5.1f}s\n{}'.format(topen, tdead, '-' * 82)
                print(line)

    def plot(self, forecast=None):
        """Plot a summary of the survey statistics.

        Requires that matplotlib is installed.
        """
        import matplotlib.pyplot as plt
        import desisurvey.plots
        assert self.validate()
        D = self._data
        S = self.scheduler
        nprograms = len(S.programs)
        npasses = len(S.pass_program)
        # Find the last day of the survey.
        last = np.argmax(np.cumsum(D['completed'].sum(axis=1))) + 1
        # Combine passes into programs.
        tsetup = np.zeros((last, nprograms))
        tsplit = np.zeros((last, nprograms))
        ntiles = np.zeros(nprograms, int)
        for passnum in range(npasses):
            pidx = S.program_index[S.pass_program[passnum]]
            tsetup[:, pidx] += D['tsetup'][:last, passnum]
            tsplit[:, pidx] += D['tsplit'][:last, passnum]
            ntiles[pidx] += np.count_nonzero(S.passnum == passnum)
        actual = np.cumsum(D['completed'], axis=0)

        dt = 1 + np.arange(len(D))
        fig, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 10))

        ax = axes[0]
        npasses = D['completed'].shape[-1]
        for pidx, program in enumerate(self.scheduler.programs):
            color = desisurvey.plots.program_color[program]
            for i, passnum in enumerate(S.program_passes[program]):
                npass = np.count_nonzero(S.passnum == passnum)
                if forecast:
                    ax.plot(dt, 100 * forecast.pass_progress[passnum] / npass, ':', c=color, lw=1)
                ax.plot(dt[:last], 100 * actual[:last, passnum] / npass,
                        lw=3, alpha=0.5, c=color, label=program if i == 0 else None)
        if forecast:
            ax.plot([], [], 'b:', lw=1, label='forecast')
        ax.legend(ncol=1)
        ax.axvline(dt[last], ls='-', c='r')
        #ax.set_xlabel('Elapsed Days')
        #ax.set_xlim(0, dt[-1] + 1)
        ax.set_ylim(0, 100)
        ax.set_ylabel('Completed [%]')
        yaxis = ax.yaxis
        yaxis.tick_right()
        yaxis.set_label_position('right')
        
        ax = axes[1]
        # Plot overheads by program.
        for pidx, program in enumerate(S.programs):
            c = desisurvey.plots.program_color[program]
            scale = 86400 / ntiles[pidx] # secs / tile
            ax.plot(dt[:last], scale * np.cumsum(tsetup[:, pidx]), '-', c=c)
            ax.plot(dt[:last], scale * np.cumsum(tsplit[:, pidx]), '--', c=c)
            ax.plot(dt[:last], scale * np.cumsum(D['tdead'][:last, pidx]), ':', c=c)
            if forecast:
                row = forecast.df.iloc[forecast.programs.index(program)]
                ax.scatter([dt[-1], dt[-1], dt[-1]], [
                    row['Setup overhead / tile (s)'],
                    row['Cosmic split overhead / tile (s)'],
                    row['Operations overhead / tile (s)']], s=50, lw=0, c=c)
        ax.plot([], [], 'b-', label='setup')
        ax.plot([], [], 'b--', label='split')
        ax.plot([], [], 'b:', label='dead')
        for pidx, program in enumerate(S.programs):
            ax.plot([], [], '-', c=desisurvey.plots.program_color[program], label=program)
        ax.legend(ncol=2)
        ax.axvline(dt[last], ls='-', c='r')
        ax.set_xlabel('Elapsed Days')
        ax.set_ylabel('Overhead / Tile [s]')
        ax.set_xlim(0, dt[-1] + 1)
        ax.set_ylim(0, None)
        yaxis = ax.yaxis
        yaxis.set_minor_locator(plt.MultipleLocator(10))
        yaxis.tick_right()
        yaxis.set_label_position('right')
        plt.subplots_adjust(hspace=0.05)
        return fig, axes