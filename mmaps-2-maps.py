#!/bin/env python
import sys
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="Convert (pseudo-)binary mmaps results to maps equivalent.")
parser.add_argument('-l', required=True, type=int, help="the column number corresponding to the concentration of the left end, the 1-x.")
parser.add_argument('-r', required=True, type=int, help="the column number corresponding to the concentration of the right end, the x.")
args = parser.parse_args()

gs = pd.read_table('gs.out', sep='\s+', header=None)
col_to_drop = list(range(gs.shape[1] - 4))
col_to_drop.remove(args.r)
multiplier = 1/(gs[args.l] + gs[args.r])[0]
gs = gs.drop(col_to_drop + [5], axis=1).sort_values(6)
gs[args.r] = gs[args.r] * multiplier
gs.to_csv('gs.out', sep=' ', header=False, index=False, float_format='%.6f')

fit = pd.read_table('fit.out', sep='\s+', header=None).drop(col_to_drop, axis=1)
fit[args.r] = fit[args.r] * multiplier
fit.to_csv('fit.out', sep=' ', header=False, index=False, float_format='%.6f')

try:
    predstr = pd.read_table('predstr.out', sep='\s+', header=None).drop(col_to_drop, axis=1)
    predstr[args.r] = predstr[args.r] * multiplier
    predstr['E'] = 0
    predstr[[args.r, 'E', 3, 4, 5]].to_csv('predstr.out', sep=' ', header=False, index=False, float_format='%.6f')
except:
    pass
