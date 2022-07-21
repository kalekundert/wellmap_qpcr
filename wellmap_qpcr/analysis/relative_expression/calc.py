#!/usr/bin/env python3

import numpy as np
import pandas as pd

"""
Typical use:

    >>> import wellmap
    >>> from wellmap_qpcr import load_cq, agg_cq, calc_Δcq
    >>> df, extras = wellmap.load(
    ...         data_loader=load_cq,
    ... )
    >>> df_cq = df.groupby(['gene', ...]).apply(agg_cq)
    >>> df_Δcq = calc_Δcq(
    ...         df_cq.loc['expt'],  # These location refer to the 'gene'
    ...         df_cq.loc['ref'],   # group, because it's first.
    ... )
"""

def agg_cq(df):
    row = pd.Series(dtype=float)
    row['n'] = len(df['cq'])
    row['n_nan'] = df['cq'].isnull().sum()
    row['cq_mean'] = df['cq'].mean()
    row['cq_median'] = df['cq'].median()
    row['cq_min'] = df['cq'].min()
    row['cq_max'] = df['cq'].max()
    row['cq_std'] = df['cq'].std()
    return row


def calc_Δcq(df_expt, df_ref, center='mean'):
    """
    The two data frames need to have the same index, so that they can be 
    subtracted from one another.  There are a few ways to do this:

    If you previously used `groupby()`, then the data frame will have a 
    hierarchical index and you can use `loc` to select the relevant 
    experimental/reference values::

        >>> df_cq = df.groupby(['gene', ...]).apply(agg_cq)
        >>> df_Δcq = calc_Δcq(
        ...         df_cq.loc['expt'],
        ...         df_cq.loc['ref'],
        ... )

    If the above doesn't apply, you can use `query()`::

        >>> df_Δcq = calc_Δcq(
        ...         df_cq.query('gene == "expt")
        ...         df_cq.query('gene == "ref")
        ... )

    If the column you're querying is an index, you need to take some extra 
    steps.  I think there are two ways to do it:

        >>> df_Δcq = calc_Δcq(
        ...         df_cq.query('gene == "expt").droplevel('gene'),
        ...         df_cq.query('gene == "ref").droplevel('gene'),
        ... )
    
    ::

        >>> df_Δcq = calc_Δcq(
        ...         df_cq.reset_index().query('gene == "expt"),
        ...         df_cq.reset_index().query('gene == "ref"),
        ... )
    """
    return _calc_delta(
            df_expt, df_ref,
            center=center, 
            prefixes=('cq', 'Δcq'),
    )

def calc_ΔΔcq(df_expt, df_ref, center='mean'):
    return _calc_delta(
            df_expt, df_ref,
            center=center, 
            prefixes=('Δcq', 'ΔΔcq'),
    )

def _calc_delta(df_expt, df_ref, center='mean', *, prefixes):
    x, x0 = df_expt, df_ref
    cq, Δcq = prefixes

    df = pd.DataFrame(index=x.index)
    df[f'{Δcq}_mean'] = x[f'{cq}_mean'] - x0[f'{cq}_mean']
    df[f'{Δcq}_median'] = x[f'{cq}_median'] - x0[f'{cq}_median']
    # https://stats.stackexchange.com/questions/112351/standard-deviation-after-subtracting-one-mean-from-another
    # https://stats.stackexchange.com/questions/25848/how-to-sum-a-standard-deviation
    df[f'{Δcq}_std'] = np.sqrt(x[f'{cq}_std']**2 + x0[f'{cq}_std']**2)

    # Assume perfect efficiency (i.e. 2).  If the reference and target genes 
    # have very different efficiencies, I might need to use the Pfaffl method.
    df['fold_change'] = 2**(-df[f'{Δcq}_{center}'])
    df['fold_change_bound'] = 2**(-df[f'{Δcq}_{center}'] + df[f'{Δcq}_std'])
    df['fold_change_err'] = df['fold_change_bound'] - df['fold_change']

    return df
