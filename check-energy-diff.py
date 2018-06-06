#!/bin/env python
import sys
import numpy as np
import pandas as pd
import atat_module
import argparse

parser = argparse.ArgumentParser(
    description="Check if the energy discrepancy is too large (measured by multiples of std), or lies too much above the hull (measured by quantile).")
parser.add_argument('-q', required=True, type=float, help="the quantile of the energy at a certain concentration, above which to remove.")
parser.add_argument('-d', required=True, type=float, help="the multiple of std, above which to remove.")
args = parser.parse_args()

quantile = float(args.q)
std_multi = float(args.d)

fit = atat_module.get_df('fit', 'fit.out')

above_limit_list = []
for c in sorted(list(set(fit['c']))):
    fit_c = fit[np.isclose(fit['c'], c)]
    fit_c_above_limit = fit_c[fit_c['E'] > fit_c['E'].quantile(quantile)]
    above_limit_list.append(fit_c_above_limit)

above_limit = pd.concat(above_limit_list)
to_be_gone = above_limit[above_limit['E_diff'] > std_multi * fit['E_diff'].std()]
print(to_be_gone)
