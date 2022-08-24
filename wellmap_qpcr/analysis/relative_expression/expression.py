#!/usr/bin/env python3

"""\
Compare relative gene expression using the ΔΔCq equation.

Usage:
    qpcr-relative-expression <toml> [-o <path> | -O] [-v]
    qpcr-relative-expression (amp|amplification) [...]
    qpcr-relative-expression melt [...]

Subcommands:
    If no subcommands are specified, a bar chart comparing the relative 
    expression of all the genes specified in the <toml> layout will be 
    displayed.  Otherwise, the analysis associated with the given subcommand 
    will be performed.  You can pass the `-h` flag to each subcommand to get 
    more information on what exactly it does.

    All of the subcommands are related in that they expect the plate layout 
    (specified by the <toml> argument) to contain the same information.
    
Arguments:
    <toml>
        A wellmap file describing the experimental layout.  See the 'Layout' 
        section below for more information on the expected contents of this 
        file.

Options:
    -o --output <path>
        Output an image of the plot to the given path, instead of launching the 
        interactive GUI.  The file type is inferred from the file extension.  
        If the path contains a percent sign (e.g. '%.svg'), it will be replaced 
        with the base name of the <toml> path.

    -O --output-default
        Output an image of the plot to the default path.  This is equivalent to 
        specifying `--output %.svg`.

    -v --verbose
        Print the raw numbers for each step of the calculation.

Layout:
    The layout of the plate should be described using the wellmap file format.  
    For a general description of this format, refer to:

        https://wellmap.readthedocs.io

    Each well should have the following attributes defined:

        label:
            A label identifying a group of wells that will be used to calculate 
            a single relative gene expression value.  Each group must contain 4 
            conditions: experimental and reference genes (as specified by the 
            "gene" attribute), and experimental and reference treatments (as 
            specified by the "treatment" attribute).  There may be multiple 
            replicates of any of these conditions.  The calculated relative 
            expression will be between the experimental and control treatments, 
            and normalized by the reference gene.

            Currently there can only be one reference gene per group, although 
            having multiple reference genes is known to improve accuracy.  If 
            you are interested in this functionality, let me know and I'll be 
            happy to implement it.

            This label will also be displayed on the resulting plots, so pick 
            something readable.  If the label contains any Python string- 
            formatting codes (e.g. "{strain}"), they will be replaced with the 
            corresponding attributes for the well in question.  This can allow 
            the label to be built from other well attributes and specified 
            succinctly for all wells at once.

        sublabel:
            Similar to label, but identifying the gene and treatment conditions 
            within each label group.  These labels are used when making plots 
            looking at the data that go into calculating the relative gene 
            expression for each label, e.g. for the "amplification" and "melt" 
            subcommands.

        housekeeping:
            A boolean indicating whether this is the housekeeping gene (True) 
            or the experimentally interesting gene (False).  This information 
            can also be specified via the `qpcr.housekeeping.*` metadata 
            options.

        treatment:
            A boolean indicating whether this is the experimental (True) or 
            reference (False) treatment condition.  This information can also 
            be specified via the `qpcr.treatment.*` metadata options.

    The following metadata can also be provided:

        qpcr.housekeeping.true
        qpcr.housekeeping.false
            Pandas query strings that can be used to find the wells 
            corresponding to the housekeeping and experimental genes, 
            respectively.

        qpcr.treatment.true
        qpcr.treatment.false
            Pandas query strings that can be used to find the wells 
            corresponding to the experimental and reference treatments, 
            respectively.

        qpcr.order
            The order in which each label should be displayed.  This should be 
            a dictionary where the keys are labels and the values are either 
            integers or lists of 2 integers.  Every value must have the same 
            form; you cannot mix integers and lists.  If you specify integers, 
            the labels will be arranged in that order.  If you specify lists, 
            the first number will be used to split the labels between multiple 
            plots and the second will be used to arrange labels within a plot.

            By default, the order will be inferred from the plate layout, with 
            the top-leftmost wells being displayed first.

        qpcr.color
            The color that should be used for each label.  This should be a 
            dictionary where the keys are labels and the values are colors.  
            The following colors can be specified, either by name or by number:

                1. blue      2. red     3. olive    4. orange
                5. purple    6. navy    7. teal

            By default, all plots will be blue.
"""

import wellmap
import sys, docopt
import matplotlib.pyplot as plt
import pandas as pd

from .calc import agg_cq, calc_Δcq, calc_ΔΔcq
from .layout import add_labels, add_ΔΔcq_flags, init_style
from wellmap_qpcr.load import load_cq
from wellmap_qpcr.utils import plot_or_save, resolve_img_path
from color_me import ucsf
from pathlib import Path

pd.options.display.width = sys.maxsize
pd.options.display.max_rows = sys.maxsize

def main():
    # It's a bit gross to be accessing `sys.argv` directly, but there's no 
    # other way to get this behavior.  See: `kalekundert/docopt-rewrite#2`
    try:
        analysis = sys.argv[1]
    except IndexError:
        analysis = None

    if analysis in ('amp', 'amplification'):
        from .amplification import main
        return main()

    if analysis == 'melt':
        from .melt import main
        return main()

    args = docopt.docopt(__doc__)
    analysis = args.get('<analysis>')
    layout_path = Path(args['<toml>'])
    img_path = resolve_img_path(
            args['--output'],
            '%_{analysis}.svg' if analysis else f'%.svg',
            args['--output-default'],
            layout_path,
    )

    df, layout, style = load(layout_path, verbose=args['--verbose'])
    style.finalize(layout)

    with plot_or_save(layout_path, img_path):
        plot_expression(df, style)

def load(layout_path, verbose=False):
    df, extra = wellmap.load(
            layout_path,
            data_loader=load_cq,
            merge_cols=True,
            path_guess='{0.stem}',
            extras=True,
    )

    add_labels(df, extra)
    add_ΔΔcq_flags(df, extra)

    def cols(*cols):
        cols = list(cols)
        if 'treatment' not in df:
            cols.remove('treatment')
        return cols

    if verbose:
        print(df[cols('well', 'housekeeping', 'treatment', 'label', 'cq')])
        print()

    layout = df

    df = df\
            .groupby(cols('housekeeping', 'treatment', 'label'))\
            .apply(agg_cq)

    if verbose:
        print(df)
        print()

    df = calc_Δcq(
            df_expt=df.loc[0],
            df_ref=df.loc[1],
    )

    if verbose:
        print(df)
        print()

    # If the user didn't specify experimental/control treatment conditions, 
    # stop here and just report ΔCq instead of ΔΔCq.

    if 'treatment' in df.index.names:
        df = calc_ΔΔcq(
                df_expt=df.loc[1],
                df_ref=df.loc[0],
        )

        if verbose:
            print(df)

    return df, layout, init_style(extra)

def plot_expression(df, style):
    n_cols, n_bars = style.shape

    fig, axes = plt.subplots(
            1, n_cols,
            squeeze=False,
            sharey=True,
    )
    axes = axes[0,:]

    ticks = {}

    for label, row in df.iterrows():
        i, j = style.indices[label]
        y = row['fold_change']
        y_err = row['fold_change_err']
        color = style.color.get(label, ucsf.blue[0])

        axes[i].plot(
                [j, j],
                [0, y],
                linewidth=5,
                solid_capstyle='butt',
                color=color,
        )
        axes[i].errorbar(
                [j], [y], [y_err],
                color=color,
        )


        ticks.setdefault(i, {})[j] = label

    for i, ax in enumerate(axes):
        ax.set_xlim(-0.5, n_bars - 0.5)
        ax.set_ylim(0)
        ax.set_xticks(
                list(ticks[i].keys()),
                labels=list(ticks[i].values()),
                rotation='vertical',
        )
        ax.axhline(1, linestyle=':', color=ucsf.light_grey[0], zorder=-1)

    axes[0].set_ylabel('gene expression')

    fig.tight_layout()








