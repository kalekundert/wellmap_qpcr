#!/usr/bin/env python3

"""\
Plot melt curves for each experimental condition.

Usage:
    qpcr-relative-expression melt <toml> [-o <path> | -O]

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
        specifying `--output %_melt.svg`.

Looking at the melt curves is a useful (but not foolproof) way to confirm that 
the PCR reactions worked cleanly.
"""

import wellmap
import docopt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .layout import add_labels, add_ΔΔcq_flags, init_style
from wellmap_qpcr.load import load_melt
from wellmap_qpcr.utils import plot_or_save, resolve_img_path
from color_me import ucsf
from more_itertools import flatten
from pathlib import Path

def main():
    args = docopt.docopt(__doc__)
    layout_path = Path(args['<toml>'])
    img_path = resolve_img_path(
            args['--output'],
            '%_melt.svg',
            args['--output-default'],
            layout_path,
    )

    df, style = load(layout_path)
    style.finalize(df)

    with plot_or_save(layout_path, img_path):
        plot_melt_groups(df, style)

def load(layout_path):
    df, extra = wellmap.load(
            layout_path,
            data_loader=load_melt,
            merge_cols=True,
            path_guess='{0.stem}',
            extras=True,
    )
    add_labels(df, extra)
    add_ΔΔcq_flags(df, extra)

    return df, init_style(extra)

def plot_melt_groups(df, style):
    n_rows, n_cols = style.shape
    fig, axes = plt.subplots(
            n_rows, n_cols,
            squeeze=False,
            sharex=True,
            sharey=True,
            figsize=(n_cols*2 + 2, n_rows*2),
    )

    for label, g in df.groupby('label'):
        ij = style.indices[label]
        labels = plot_melt_curves(axes[ij], label, g, style)

    for ax in axes[:,0]:
        ax.set_ylabel('dRFU/dT')
    for ax in axes[-1,:]:
        ax.set_xlabel('cycles')

    axes[0,-1].legend(
            labels.values(), labels.keys(),
            bbox_to_anchor=(1,1),
            loc='upper left',
    )

    fig.tight_layout()

def plot_melt_curves(ax, label, df, style):
    ax.set_title(label)

    labels = {}
    cols = ['well', 'sublabel', 'gene', 'treatment']

    for (well, sublabel, gene, treatment), g in df.groupby(cols):
        color = style.color.get(sublabel, ucsf.blue[0]) \
                if treatment else ucsf.dark_grey[0]
        linestyle = '-' if gene else '--'

        artists = ax.plot(
                g['temp_C'], g['rfu_deriv'],
                label=sublabel,
                color=color,
                linestyle=linestyle,
                zorder=treatment,
        )
        labels[sublabel] = artists[0]

    return labels

