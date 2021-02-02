#!/usr/bin/env python3

#from slugify import slugify
import pandas as pd
from functools import partial

def load_cq(path):
    if path.is_dir():
        path = path / 'Quantification Cq Results.csv'

    loaders = {
            '.csv': pd.read_csv,
            '.tsv': partial(pd.read_csv, sep='\t'),
            '.xslx': pd.read_excel,
    }
    return loaders[path.suffix](path)\
            .rename(columns={'Well': 'well0', 'Cq': 'cq'})

def compact(dir):
    pass


def _find_biorad_files(dir):
    for src in dir.glob('* -  *'):
        m = re.match(r'(src.name.')




