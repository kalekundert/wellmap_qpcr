#!/usr/bin/env python3

import pandas as pd
from functools import partial
from more_itertools import one
from os.path import getmtime

LOADERS = {
        '.csv': pd.read_csv,
        '.tsv': partial(pd.read_csv, sep='\t'),
        '.xslx': pd.read_excel,
}

def load_cq(path):
    if path.is_dir():
        path = one(path.glob('Quantification Cq Results.*'))

    return LOADERS[path.suffix](path)\
            .rename(columns={'Well': 'well0', 'Cq': 'cq'})

def load_trace(path):
    if path.is_dir():
        path = one(path.glob('Quantification Amplification Results*'))

    return LOADERS[path.suffix](path)\
            .rename(columns={'Cycle': 'cycle'})\
            .melt(id_vars=['cycle'], var_name='well', value_name='rfu')

def load_melt(path):
    if path.is_dir():
        path = one(path.glob('Melt Curve Derivative Results*'))

    return LOADERS[path.suffix](path)\
            .rename(columns={'Temperature': 'temp_C'})\
            .melt(id_vars=['temp_C'], var_name='well', value_name='rfu_deriv')


