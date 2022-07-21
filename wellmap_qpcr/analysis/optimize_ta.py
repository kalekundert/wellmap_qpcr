#!/usr/bin/env python3

import wellmap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from color_me.ucsf import iter_colors
from wellmap_qpcr import load_cq
from .main import App

class OptimizeTa(App):
    """\
Plot amplification for different annealing temperatures.

Usage:
    qpcr-optimize-ta <toml> [-o <path>]

Arguments:
    <toml>
        A wellmap file describing the experimental layout.  The layout should 
        contain the following information:

        For each well:

          template:
              The name of the template being amplified.  Wells with different 
              templates will be plotted separately.  If not specified, it will 
              be assumed that wells have the same template.

          primers:
              The name of the primers being used.  Wells with different primers 
              will be plotted separately.  If not specified, it will be assumed 
              that all wells use the same primers.

          anneal_temp_C:
              The annealing temperature in effect for this well, in °C.  
              Typically, you would setup a temperature gradient such that each 
              row or column has a distinct temperature.

Options:
    -o --output <path>
        Output an image of the plot to the given path, instead of launching the 
        interactive GUI.  The file type is inferred from the file extension.  
        If the path contains a dollar sign (e.g. '$.svg'), it will be replaced 
        with the base name of the <toml> path.
"""

    def __init__(self, layout_toml):
        self.layout_toml = layout_toml

    def load(self):
        df = wellmap.load(
                self.layout_toml,
                data_loader=load_cq,
                merge_cols=True,
                path_guess='{0.stem}',
        )

        # Fill in optional columns:
        if 'template' not in df:
            df['template'] = np.nan
        if 'primers' not in df:
            df['primers'] = np.nan

        return df, {}

    def plot(self, df):
        fig, ax = plt.subplots()
        groups = df.groupby(['template', 'primers'], dropna=False)
        have_labels = False

        for color, (key, g) in iter_colors(groups):
            label = ', '.join((x for x in key if not np.isnan(x)))
            have_labels = have_labels or bool(label)

            ax.plot(
                    g['anneal_temp_C'],
                    g['cq'],
                    marker='+',
                    linestyle='none',
                    color=color,
                    label=label,
            )
        
        if have_labels:
            ax.legend()

        ax.set_xlabel('$T_A$ (°C)')
        ax.set_ylabel('$C_q$', rotation='horizontal', ha='right')

        fig.tight_layout()

        return fig
