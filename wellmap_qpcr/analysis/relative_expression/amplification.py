#!/usr/bin/env python3

"""\
Plot raw amplification curves for each experimental condition.

This can reveal information that would be abstracted away in the relative 
expression bar plot, e.g. whether the amplification was clean, whether certain 
conditions have abnormally high/low Cq values, etc.

Usage:
    qpcr-relative-expression (amp|amplification) <toml> [-o <path> | -O] [-l]

Arguments:
    <toml>
        A wellmap file describing the experimental layout.  Refer to 
        `qpcr-relative-expression -h` for a detailed description of this file.

Options:
    -o --output <path>
        Output an image of the plot to the given path, instead of launching the 
        interactive GUI.  The file type is inferred from the file extension.  
        If the path contains a percent sign (e.g. '%.svg'), it will be replaced 
        with the base name of the <toml> path.

    -O --output-default
        Output an image of the plot to the default path.  This is equivalent to 
        specifying `--output %_amp.svg`.

    -l --log-rfu
        Plot the relative fluorescence unit (RFU) axis on a log scale.
"""

import wellmap
import docopt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .layout import add_labels, add_ΔΔcq_flags, init_style
from wellmap_qpcr.load import load_trace, load_cq
from wellmap_qpcr.utils import plot_or_save, resolve_img_path
from color_me import ucsf
from more_itertools import flatten
from pathlib import Path

def main():
    args = docopt.docopt(__doc__)
    layout_path = Path(args['<toml>'])
    img_path = resolve_img_path(
            args['--output'],
            '%_amp.svg',
            args['--output-default'],
            layout_path,
    )

    df_cq, df_trace, style = load(layout_path)
    style.finalize(df_cq)
    
    with plot_or_save(layout_path, img_path):
        plot_trace_groups(df_cq, df_trace, style, args['--log-rfu'])

def load(layout_path):
    layout, extra = wellmap.load(
            layout_path,
            path_guess='{0.stem}',
            path_required=True,
            extras=True,
    )

    add_labels(layout, extra)
    add_ΔΔcq_flags(layout, extra)

    # The `load_data()` function should probably be provided by wellmap...
    cq = load_data(layout, load_cq)
    trace = load_data(layout, load_trace)

    # Wellmap should probably also provided a function that does the merge with 
    # the same API/semantics as `wellmap.load()`...
    df_cq = pd.merge(layout, cq, on=['well0'])
    df_trace = pd.merge(layout, trace, on=['well'])

    return df_cq, df_trace, init_style(extra)

def load_data(layout, data_loader):
    chunks = []

    for path in layout['path'].unique():
        df = data_loader(path)
        df['path'] = path
        chunks.append(df)

    return pd.concat(chunks, sort=False)

def plot_trace_groups(df_cq, df_trace, style, log_rfu=False):
    n_rows, n_cols = style.shape
    fig, axes = plt.subplots(
            n_rows, n_cols,
            squeeze=False,
            sharex=True,
            sharey=True,
            figsize=(n_cols*2 + 2, n_rows*2),
    )

    for label, g in df_trace.groupby('label'):
        ij = style.indices[label]
        labels = plot_trace_group(
                axes[ij], label,
                df_cq.query('label==@label'), g, 
                style,
        )

    for ax in axes[:,0]:
        ax.set_ylabel('RFU')
    for ax in axes[-1,:]:
        ax.set_xlabel('cycles')

    for ax in axes.flat:
        ax.set_yscale('symlog' if log_rfu else 'linear')

    axes[0,-1].legend(
            labels.values(), labels.keys(),
            bbox_to_anchor=(1,1),
            loc='upper left',
    )

    fig.tight_layout()

def plot_trace_group(ax, label, df_cq, df_trace, style):
    ax.set_title(label)

    cols = ['well', 'sublabel', 'housekeeping', 'treatment']
    df_cq = df_cq.set_index('well')

    if 'treatment' not in df_trace:
        df_trace['treatment'] = True

    labels = {}

    for (well, sublabel, housekeeping, treatment), g in df_trace.groupby(cols):
        color = style.color.get(sublabel, ucsf.blue[0]) \
                if treatment else ucsf.dark_grey[0]
        linestyle = '--' if housekeeping else '-'
        marker = '+' if housekeeping else 'o'

        artists = ax.plot(
                g['cycle'], g['rfu'],
                label=sublabel,
                color=color,
                linestyle=linestyle,
                zorder=treatment,
        )
        labels[sublabel] = artists[0]

        cq = df_cq.loc[well,'cq']
        rfu = np.interp(cq, g['cycle'], g['rfu'])
        ax.plot(
                [cq], [rfu], 
                marker=marker,
                markeredgecolor=color,
                markerfacecolor='none',
                zorder=treatment+2,
        )

    return labels

