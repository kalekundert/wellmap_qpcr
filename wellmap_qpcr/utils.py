#!/usr/bin/env python3

import wellmap
import matplotlib.pyplot as plt
import os, sys

from pathlib import Path
from contextlib import contextmanager

# Way to plot or savefig.  Need to os.fork before plotting.  Want to use with 
# statement, save gcf().  Good enough for now, I can be more clever later.

@contextmanager
def plot_or_save(layout_path, img_path, fork=True):
    if fork and not img_path:
        if os.fork() != 0:
            sys.exit()

    yield

    if img_path:
        plt.savefig(img_path)
        plt.close()
    else:
        cmd = Path(sys.argv[0]).name
        plt.gcf().canvas.set_window_title(f'{cmd} {layout_path}')
        plt.show()

def resolve_img_path(img_template, default_img_template, use_default, layout_path):
    if not img_template and not use_default:
        return None
    if use_default:
        img_template = default_img_template

    return Path(img_template.replace('%', layout_path.stem))

def parse_wells(well_strs):
    ijs = flatten([
        wellmap.iter_well_indices(x)
        for x in well_strs
    ])
    return [
            wellmap.well_from_ij(*ij)
            for ij in ijs
    ]

def filter_by_label(labels, df):
    pass

def filter_by_wells(wells, df):
    return df[df['well'].isin(wells)]

