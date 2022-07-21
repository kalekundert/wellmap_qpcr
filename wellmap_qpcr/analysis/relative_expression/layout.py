#!/usr/bin/env python3

import numpy as np

from pydantic import BaseModel, validator
from typing import Optional, Dict, Tuple
from more_itertools import unique_everseen as unique
from color_me import ucsf

# The style class is nice for the common [extra] fields (e.g.  order, color, 
# etc.), but how to get custom fields?

class Style(BaseModel):
    order: Dict[str,int] | Dict[str,Tuple[int,int]] = {}
    color: Dict[str,str] = {}
    shape: Optional[Tuple[int,int]] = None
    indices: Optional[Dict[str,Tuple[int,int]]] = None

    @validator('color', pre=True)
    def parse_color(cls, colors):
        return {k: parse_color(v) for k, v in colors.items()}

    def finalize(self, df):
        """
        Update the style to reflect the layout being plotted.
        """

        if not self.order:
            self.order = infer_order_from_layout(df)

        if all(isinstance(x, int) for x in self.order):
            self.shape, self.indices = infer_shape_from_order(self.order)
        else:
            self.shape, self.indices = infer_shape_from_order_tuples(self.order)

def add_labels(df, extra):
    # - "label" is used to identify a group of wells that can be collectively 
    #   used to calculate a ΔΔCq value, i.e. relative gene expression.  Each 
    #   group of wells with the same label must have at least 4 conditions: 
    #   experimental/ref gene/treatment.
    # - "sublabel" is used to identify the wells within a "label" group, for 
    #   cases where the "label" group is clear from context (e.g. each group is 
    #   split into its own plot).

    format_cols = ['label', 'sublabel']

    def format_cell(row, col):
        return row[col].format_map(row)
    
    for col in format_cols:
        df[col] = df.apply(format_cell, args=(col,), axis=1)

def add_ΔΔcq_flags(df, extra):
    for col in ['gene', 'treatment']:
        if col in df.columns:
            continue

        i = df.eval(extra['qpcr'][col]['expt'])
        j = df.eval(extra['qpcr'][col]['ref'])

        df[col] = np.nan
        df.loc[i,col] = True
        df.loc[j,col] = False

def init_style(extra):
    return Style.parse_obj(extra.get('qpcr', {}))

def infer_order_from_layout(df):
    order = {}
    if 'row_i' in df and 'col_j' in df:
        for label, g in df.groupby('label'):
            order[label] = min(
                    (x['row_i'], x['col_j'])
                    for i, x in g.iterrows()
            )
    return order

def infer_shape_from_order(order):
    values = sorted(unique(order))
    shape = 0, len(values)

    index_map = {k: v for v, k in enumerate(values)}
    indices = {k: index_map[v] for k, v in order.items()}

    return shape, indices

def infer_shape_from_order_tuples(order):
    xs, ys = zip(*order.values())
    xs = sorted(unique(xs))
    ys = sorted(unique(ys))

    shape = len(xs), len(ys)

    x_index_map = {k: v for v, k in enumerate(xs)}
    y_index_map = {k: v for v, k in enumerate(ys)}
    indices = {
            k: (x_index_map[x], y_index_map[y])
            for k, (x, y) in order.items()
    }

    return shape, indices

def parse_color(color):
    if isinstance(color, int):
        return ucsf.cycle[color-1]
    else:
        return getattr(ucsf, color)


