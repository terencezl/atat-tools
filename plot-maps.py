#!/bin/env python
from scipy import interpolate
import pandas as pd
import argparse
import atat_module
import matplotlib as mpl

parser = argparse.ArgumentParser(
    description="Plot the ground state analysis result with maps convex hull.")
parser.add_argument('--gs', action='store_true', help='include the gs analysis structures.')
parser.add_argument('--gsline-vasp', action='store_true', help='use the VASP calculated values to form the gs line.')
parser.add_argument('--gslim', nargs=2, metavar=('from', 'to'), type=float, help="the range of ground states, from 0 to 1. The ref energy will be adjusted")
parser.add_argument('--xlim', nargs=2, metavar=('from', 'to'), type=float, help="the range of x-axis")
parser.add_argument('--ylim', nargs=2, metavar=('from', 'to'), type=float, help="the range of y-axis")
parser.add_argument('-s', metavar='savefig_name', help="if specified, save to file but not display")
args = parser.parse_args()

if args.s:
    mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('research')

#plt.rcParams['axes.prop_cycle'] = plt.rcParamsDefault['axes.prop_cycle']
colors = plt.rcParams['axes.color_cycle']

if args.gs:
    strcorrenergy = pd.read_csv('strcorrenergy.csv', index_col='str_id')
predstr = atat_module.get_df('predstr', 'predstr.out')
fit = atat_module.get_df('fit', 'fit.out')
gs = atat_module.get_df('gs', 'gs.out')

if args.gsline_vasp:
    gsline_E = '-'
    gsline_E_fit = ''
    interp_func = interpolate.interp1d(gs['c'], gs['E'])
else:
    gsline_E = ''
    gsline_E_fit = '-'
    interp_func = interpolate.interp1d(gs['c'], gs['E_fit'])

if args.gslim:
    x1 = args.gslim[0]
    x2 = args.gslim[1]
    # gs = gs[(gs['c'] >= args.gslim[0]) & (gs['c'] <= args.gslim[1])]

    def E_ref(x):
        return (x2 - x)/(x2 - x1) * interp_func(x1) + (x - x1)/(x2 - x1) * interp_func(x2)

    if args.gs:
        strcorrenergy['E_per_site'] = strcorrenergy['E_per_site'] - E_ref(strcorrenergy['c'])
    predstr['E_fit'] = predstr['E_fit'] - E_ref(predstr['c'])
    fit['E_fit'] = fit['E_fit'] - E_ref(fit['c'])
    fit['E'] = fit['E'] - E_ref(fit['c'])
    gs['E_fit'] = gs['E_fit'] - E_ref(gs['c'])
    gs['E'] = gs['E'] - E_ref(gs['c'])

# plot
if args.gs:
    plt.plot(strcorrenergy['c'], strcorrenergy['E_per_site'], '+', color=colors[2], label='genstr.out', alpha=0.7)
plt.plot(predstr['c'], predstr['E_fit'], 'x', color=colors[1], label='predstr.out', alpha=0.7)
plt.plot(fit['c'], fit['E_fit'], '+', ms=7, color=colors[0], mew=1.2, label='fit.out', alpha=1)
plt.plot(fit['c'], fit['E'], 'o', ms=7, mfc='none', mec=colors[0], mew=1.2, label='fit.out', alpha=1)
plt.plot(gs['c'], gs['E_fit'], gsline_E_fit + '+', ms=7, color='k', mew=1.2, label='gs.out', alpha=1)
plt.plot(gs['c'], gs['E'], gsline_E + 'o', ms=7, color='k', mfc='none', mec='k', mew=1.2, label='gs.out', alpha=1)
plt.xlabel('Concentration ($x$)')
plt.ylabel('Energy (eV)')
ax = plt.gca()
if args.gs:
    plt.legend([(ax.lines[4], ax.lines[5]), (ax.lines[2], ax.lines[3]), ax.lines[1], ax.lines[0]], ['ground states', 'known structures', 'predicted structures', 'enumerated structures'], fontsize='x-small', frameon=True, framealpha=0.7, loc=0)
else:
    plt.legend([(ax.lines[3], ax.lines[4]), (ax.lines[1], ax.lines[2]), ax.lines[0]], ['ground states', 'known structures', 'predicted structures'], fontsize='small', frameon=True, framealpha=0.7, loc=0)

if args.xlim:
    plt.xlim(*args.xlim)
if args.ylim:
    plt.ylim(*args.ylim)

plt.tight_layout()

if args.s:
    plt.savefig(args.s)
else:
    plt.show()
