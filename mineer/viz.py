"""
Functions for displaying visulizations
"""

import functools
import os
from .utils import Project
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_context('notebook', font_scale=1.25)

class Viz:
    def __init__(self, project: Project, outdir: str, nreads=5000):
        """
        Project visualizer
        nreads: number of reads to subsample from existing subsample for ploting
        """
        self.project = project
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        self.outdir = outdir
        if nreads > project.nreads:
            nreads = project.nreads
        self.nreads = nreads
    
        self.kwargs = dict(col='direction', col_order=['f', 'r'], height=5, aspect=1.5)
    @functools.cached_property
    def subreads(self):
        srs = self.project.fwd_sub[:self.nreads]
        if self.project.paired:
            srs.extend(self.project.rev_sub[:self.nreads])
        return srs

    @property
    def untrimmed_phred_df(self):
        """DataFrame of PHRED at each position of subsampled reads"""
        max_len = max([len(r.untrimmed.record) for r in self.subreads])
        nsubreads = len(self.subreads)
        phreds = np.empty((nsubreads, max_len))
        phreds[:] = np.nan
        directions = [None] * nsubreads

        for i in range(nsubreads):
            r = self.subreads[i]
            # Add direction
            directions[i] = r.file.direction
            # Add phred
            phred = r.untrimmed.phred
            length = r.untrimmed.length
            phreds[i, :length] = phred
        phred_df = pd.DataFrame(phreds)
        # Add direction
        phred_df = phred_df.join(pd.Series(directions).rename('direction')).rename_axis(index='read', columns='pos')
        
        return phred_df

    @property
    def plotPhredProfiles(self):
        # Wideframe dataframe of PHRED scores for each subsampled read with directions 
        phred_df = self.untrimmed_phred_df
        # Melt for plotting boxplots
        melted = phred_df.reset_index().melt(['read', 'direction'], value_name='phred')
        max_len = melted.pos.max()

        g = sns.FacetGrid(melted, **self.kwargs)
        g.map(sns.boxplot, 'pos', 'phred', color='.1', showfliers=False, order=range(max_len), width=.1, whis=1)
        for ax, direction in zip(g.axes.ravel(), g.col_names):
            data = phred_df[phred_df.direction.eq(direction)]
            ax.plot(data.mean(numeric_only=True), lw=1.5, label='Mean')
            ax.plot(data.median(numeric_only=True), lw=1.5, label='Median')
        ax.legend(bbox_to_anchor=(1, 1))
        g.set(xticks=range(max_len)[::20])
        g.set_xlabels('Position')
        g.set_ylabels('PHRED score')
        return g

    @property
    def plotTruncPosDist(self):
        trimpos_df = pd.DataFrame([{'trimstart': r.trimpos_mineer[0], 'trimend': r.trimpos_mineer[1], 'direction': r.file.direction} for r in self.subreads])
        trimpos_melted = trimpos_df.melt('direction', var_name='position', value_name='pos')
        # Plot distributions
        g = sns.displot(trimpos_melted, x='pos', hue='position', palette='Set1', **self.kwargs)
        # Add lines for global truncation positions
        for ax, trimpos in zip(g.axes.ravel(), [self.project.fwd_pos, self.project.rev_pos]):
            for p in trimpos:
                ax.axvline(p, color='gold', lw=4)

        g.set_xlabels('Position')
        return g

    def genFigs(self):
        # Generate phred plot
        phred_filepath = os.path.join(self.outdir, 'phred_profiles')
        self.plotPhredProfiles
        plt.savefig(phred_filepath, facecolor='w', bbox_inches='tight')
        plt.close()

        # Generate truncation position plot
        trun_filepath = os.path.join(self.outdir, 'trunc_dist')
        self.plotTruncPosDist
        plt.savefig(trun_filepath, facecolor='w', bbox_inches='tight')
        plt.close()