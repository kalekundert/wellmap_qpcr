#!/usr/bin/env python3

import wellmap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import autoprop

from color_me import ucsf
from scipy.stats import linregress
from wellmap_qpcr import load_cq
from matplotlib.lines import Line2D
from dataclasses import dataclass, fields
from more_itertools import mark_ends
from datetime import datetime
from numpy import log10
from .main import App

@autoprop
class CheckEfficiency(App):
    """\
Plot and analyze qPCR standard curves.

Usage:
    qpcr-check-efficiency <toml> [-o <path>]

Arguments:
    <toml>
        A wellmap file describing the experimental layout.  The layout should 
        contain the following parameters (<well> indicates any well 
        specification, e.g. `row.A`, `col.1`, `block.2x2.A1`, etc.):

        <well>.template_conc: number, required
            The concentration of the template in the given well.  Any well 
            missing this information will not be included in the plot or the 
            analysis.  See also: `qpcr.conc_unit`.

        <well>.control: string, optional
            The name of the control in the given well, e.g. "No RT".  Control 
            wells are excluded from the efficiency analysis, and are labeled 
            separately on the plot.

        <well>.template: string, optional
            The name of the template used in the given well.  Wells with 
            different templates will be plotted and analyzed separately.

        <well>.primers: string, optional
            The name of the primers used in the given well.  Wells with 
            different primers will be plotted separately.

        <well>.date: date, optional
            The date of the experiment.  Experiments done on different dates 
            will be plotted separately.
              
        qpcr.conc_unit: string, optional
            The unit for the concentrations given for each well.  If given, 
            this information is added to the axis label.

Options:
    -o --output <path>
        Output an image of the plot to the given path, instead of launching the 
        interactive GUI.  The file type is inferred from the file extension.  
        If the path contains a dollar sign (e.g. '$.svg'), it will be replaced 
        with the base name of the <toml> path.

Performing a standard curve is one step in the process of validating a new pair 
of qPCR primers.  The rule of thumb is to find primers that have R²>0.99 and 
95–105% efficiency.  That said, it can be possible to account for poor primer 
efficiency in downstream analysis steps.
"""
    
    def __bareinit__(self):
        self._df = None
        self._extras = None

    def plot(self, fig_factory=plt.subplots):
        df, extras = self.df, self.extras
        expts = ['template', 'primers', 'date']
        expt_groups = [
                (ExptKey(*expt), g)
                for expt, g in df.groupby(expts, dropna=False)
        ]
        expt_colors = {
                expt: color[0]
                for color, (expt, g) in ucsf.iter_colors(expt_groups)
        }

        control_df = df[df['is_control']]
        control_groups = [
                (control, ExptKey(*expt), g)
                for (control, *expt), g in control_df.groupby(
                    ['control', *expts],
                    dropna=False,
                )
        ]
        control_markers = {
                k: (3 if i == 0 else i+4, 2, 0)
                for i, (k, g) in enumerate(control_df.groupby(['control']))
        }

        fig, ax = fig_factory()
        ax.set_title(self.layout_toml)
        legend_artists = {}

        x = df['template_conc']
        x_lim = Limits(min(x[x > 0]), max(x))

        # Plot the standard curves:

        for expt, g in expt_groups:
            i = ~g['is_control']
            x, y = g['template_conc'][i], g['cq'][i]; j = x > 0
            m, b, r, p, err = linregress(log10(x[j]), y[j])

            x_fit = np.logspace(log10(x_lim.min), log10(2*x_lim.max))
            y_fit = np.polyval((m, b), log10(x_fit))

            # The efficiency calculation will be a little different if the y 
            # values are dilutions instead of concentrations.  Consult the 
            # derivation in `docs/efficiency.lyx`.
            eff = 100 * (10**(-1/m) - 1)

            color = expt_colors[expt]
            marker = (4, 2, 0)
            label = '\n'.join([
                *expt.labels,
                f'R²={r**2:.5f}',
                f'eff={eff:.2f}%',
            ])

            ax.plot(
                    x, y,
                    color=color,
                    marker=marker,
                    linestyle='none',
            )
            ax.plot(
                    x_fit, y_fit, '--',
                    color=color,
            )

            legend_artists[expt] = [
                    Line2D(
                        [], [],
                        color=color,
                        marker=marker,
                        linestyle='--',
                        label=label,
                    ),
            ]

        # Plot the controls:

        for control, expt, g in control_groups:
            artist, = ax.plot(
                    g['template_conc'],
                    g['cq'],
                    color=expt_colors[expt],
                    marker=control_markers[control],
                    linestyle='none',
                    label=control,
            )
            legend_artists[expt].append(artist)

        x_label = '[template]'
        if 'conc_unit' in extras:
            x_label += f' ({extras["conc_unit"]})'

        ax.set_xscale(
                'symlog',
                linthresh=x_lim.min,
                subs=range(2, 10),
        )
        ax.set_xlim(-x_lim.min / 5, 2 * x_lim.max)
        ax.set_ylim(0, 40)
        ax.set_xlabel(x_label)
        ax.set_ylabel('$C_q$', rotation='horizontal', ha='right')

        handles = []
        for is_first, is_last, artists in mark_ends(legend_artists.values()):
            handles += artists

            # Add some blank space between labels.
            if not is_last:
                handles.append(Line2D([], [], linestyle='none'))

        legend = ax.legend(
                handles=handles,
                bbox_to_anchor=(1.04, 1.00),
                loc='upper left',
                borderaxespad=0,
        )
        fig.tight_layout()

        return fig

    def get_df(self):
        if self._df is None:
            self._load()
        return self._df

    def set_df(self, df):
        self._df = df

    def get_extras(self):
        if self._extras is None:
            self._load()
        return self._extras

    def set_extras(self, extras):
        self._extras = extras

    def _load(self):
        df, extras = wellmap.load(
                self.layout_toml,
                data_loader=load_cq,
                merge_cols=True,
                path_guess='{0.stem}',
                extras=['qpcr'],
        )

        def fill_default(k, default=pd.NA):
            df[k] = df[k].fillna(pd.NA) if k in df else default

        fill_default('control')
        fill_default('template')
        fill_default('primers')
        fill_default('date')

        df['is_control'] = df['control'].fillna(False).astype(bool)

        self._df = df
        self._extras = extras



@dataclass(eq=True, frozen=True)
class ExptKey:
    template: str
    primers: str
    date: datetime
    
    def __post_init__(self):
        # Replace missing values with None, so that equality comparisons will 
        # work as expected.
        for field in fields(self):
            if pd.isnull(getattr(self, field.name)):
                object.__setattr__(self, field.name, None)

    @property
    def labels(self):
        return filter(bool, [
                self.template,
                self.primers,
                self.date.strftime('%b %-d, %Y') if self.date else None,
        ])

@dataclass
class Limits:
    min: float
    max: float

### TODO: convert this script to a functional style.

def main():
    args = docopt.docopt(__doc__)
    layout_path = Path(args['<toml>'])
    img_path = resolve_img_path(
            args['--output'],
            '%.svg',
            args['--output-default'],
            layout_path,
    )

    df, style = load(layout_path)
    style.finalize(df)

    with plot_or_save(layout_path, img_path):
        plot_efficiency(df, style)

def load(layout_path):
    df, extras = wellmap.load(
            layout_path,
            data_loader=load_cq,
            merge_cols=True,
            path_guess='{0.stem}',
            extras=True,
    )

    # Should group by label instead of template/primers/date/etc.

    def fill_default(k, default=pd.NA):
        df[k] = df[k].fillna(pd.NA) if k in df else default

    fill_default('control')
    fill_default('template')
    fill_default('primers')
    fill_default('date')

    df['is_control'] = df['control'].fillna(False).astype(bool)

    return df, init_style(extra)

def plot_efficiency(df, style):
    expt_groups = df.groupby('label', dropna=False)
    expt_colors = {
            expt: color[0]
            for color, (expt, g) in ucsf.iter_colors(expt_groups)
    }

    control_df = df[df['is_control']]
    control_groups = [
            (control, label, g)
            for (control, *expt), g in control_df.groupby(
                ['control', *expts],
                dropna=False,
            )
    ]
    control_markers = {
            k: (3 if i == 0 else i+4, 2, 0)
            for i, (k, g) in enumerate(control_df.groupby(['control']))
    }

    fig, ax = plt.subplots()
    legend_artists = {}

    x = df['template_conc']
    x_lim = Limits(min(x[x > 0]), max(x))

    # Plot the standard curves:

    for expt, g in expt_groups:
        i = ~g['is_control']
        x, y = g['template_conc'][i], g['cq'][i]; j = x > 0
        m, b, r, p, err = linregress(log10(x[j]), y[j])

        x_fit = np.logspace(log10(x_lim.min), log10(2*x_lim.max))
        y_fit = np.polyval((m, b), log10(x_fit))

        # The efficiency calculation will be a little different if the y 
        # values are dilutions instead of concentrations.  Consult the 
        # derivation in `docs/efficiency.lyx`.
        eff = 100 * (10**(-1/m) - 1)

        color = expt_colors[expt]
        marker = (4, 2, 0)
        label = '\n'.join([
            *expt.labels,
            f'R²={r**2:.5f}',
            f'eff={eff:.2f}%',
        ])

        ax.plot(
                x, y,
                color=color,
                marker=marker,
                linestyle='none',
        )
        ax.plot(
                x_fit, y_fit, '--',
                color=color,
        )

        legend_artists[expt] = [
                Line2D(
                    [], [],
                    color=color,
                    marker=marker,
                    linestyle='--',
                    label=label,
                ),
        ]

    # Plot the controls:

    for control, expt, g in control_groups:
        artist, = ax.plot(
                g['template_conc'],
                g['cq'],
                color=expt_colors[expt],
                marker=control_markers[control],
                linestyle='none',
                label=control,
        )
        legend_artists[expt].append(artist)

    x_label = '[template]'
    if 'conc_unit' in extras:
        x_label += f' ({extras["conc_unit"]})'

    ax.set_xscale(
            'symlog',
            linthresh=x_lim.min,
            subs=range(2, 10),
    )
    ax.set_xlim(-x_lim.min / 5, 2 * x_lim.max)
    ax.set_ylim(0, 40)
    ax.set_xlabel(x_label)
    ax.set_ylabel('$C_q$', rotation='horizontal', ha='right')

    handles = []
    for _, is_last, artists in mark_ends(legend_artists.values()):
        handles += artists

        # Add some blank space between labels.
        if not is_last:
            handles.append(Line2D([], [], linestyle='none'))

    legend = ax.legend(
            handles=handles,
            bbox_to_anchor=(1.04, 1.00),
            loc='upper left',
            borderaxespad=0,
    )
    fig.tight_layout()

    return fig

