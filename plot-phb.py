#!/bin/env python

import argparse
import atat_module
import matplotlib as mpl

parser = argparse.ArgumentParser(
    description="Plot the T-x curves in the format of phb output.")
parser.add_argument('file', nargs='+', help="plot from FILEs (should be 1 or more)")
parser.add_argument('--xlim', nargs=2, metavar=('from', 'to'), type=float, help="the range of x-axis")
parser.add_argument('--ylim', nargs=2, metavar=('from', 'to'), type=float, help="the range of y-axis")
parser.add_argument('-s', metavar='savefig_name', help="if specified, save to file but not display")
args = parser.parse_args()

if args.s:
    mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('research')

colors = plt.rcParams['axes.color_cycle']
n_colors = len(colors)

for idx, i in enumerate(args.file):
    df = atat_module.get_df('phb', i)
    color = colors[idx % n_colors]
    plt.plot(df['x1'], df.index.values, 'x-', ms=3, color=color, label=i)
    plt.plot(df['x2'], df.index.values, 'x-', ms=3, color=color)

plt.xlabel(r'Concentration (-1 ~ 1)')
plt.ylabel('Temperature (K)')
if args.xlim:
    plt.xlim(*args.xlim)
if args.ylim:
    plt.ylim(*args.ylim)
plt.legend(loc=0, fontsize='x-small', handlelength=1, labelspacing=0.1)
plt.tight_layout()

if args.s:
    plt.savefig(args.s)
else:
    plt.show()
