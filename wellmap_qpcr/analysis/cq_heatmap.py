#!/usr/bin/env python3

import wellmap
import matplotlib.pyplot as plt
import autoprop

from .main import App
from ..load import load_cq
from wellmap import row_from_i, col_from_j, inclusive_range

@autoprop
class CqHeatmap(App):
    """\
Plot the Cq value of each reaction in the given experiment.

Usage:
    qpcr-cq-heatmap <toml> [-o <path>]

Arguments:
    <toml>
        A wellmap file describing the experimental layout.

Options:
    -o --output <path>
        Output an image of the plot to the given path, instead of launching the 
        interactive GUI.  The file type is inferred from the file extension.  
        If the path contains a dollar sign (e.g. '$.svg'), it will be replaced 
        with the base name of the <toml> path.
"""

    def __bareinit__(self):
        self._df = None

    def plot(self, fig_factory=plt.subplots):
        df = self.df
        fig, ax = fig_factory()
        img = df.pivot(index='row', columns='col', values='cq')

        artist = ax.imshow(img.values)
        plt.colorbar(artist, ax=ax, label='Cq')

        row_0 = df['row_i'].min()
        row_n = df['row_i'].max()
        num_rows = row_n - row_0 + 1
        rows = map(row_from_i, inclusive_range(row_0, row_n))

        col_0 = df['col_j'].min()
        col_n = df['col_j'].max()
        num_cols = col_n - col_0 + 1
        cols = map(col_from_j, inclusive_range(col_0, col_n))

        ax.set_yticks(range(num_rows), rows)
        ax.set_xticks(range(num_cols), cols)

        fig.tight_layout()

        return fig

    def get_df(self):
        if self._df is None:
            self._df = wellmap.load(
                    self.layout_toml,
                    data_loader=load_cq,
                    merge_cols=True,
                    path_guess='{0.stem}',
            )
        return self._df

    def set_df(self, df):
        self._df = df

