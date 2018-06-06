#!/bin/env python

from subprocess import check_output
from io import StringIO
import argparse
import numpy as np
import atat_module
import matplotlib as mpl

parser = argparse.ArgumentParser(
    description="Plot the cluster ECIs in the format of clusinfo.out.")
parser.add_argument('file', nargs='?', help="if specified, plot from file")
parser.add_argument('-s', metavar='savefig_name', help="if specified, save to file but not display")
args = parser.parse_args()

if args.s:
    mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('research')

if args.file:
    f = args.file
else:
    f = StringIO(check_output("getclus -e | grep -v '^0' | grep -v '^1'", shell=True).decode())

clusinfo = atat_module.get_df('clusinfo', f)
clusinfo['eci'] *= 1e3
pair = clusinfo[clusinfo['n_pts'] == 2]
trip = clusinfo[clusinfo['n_pts'] == 3]
quad = clusinfo[clusinfo['n_pts'] == 4]

is_trip = None
is_quad = None

if trip.values.any():
    is_trip = True
    orig_for_trip = np.ceil(pair['d'].iloc[-1]) + 1
if quad.values.any():
    is_quad = True
    orig_for_quad = orig_for_trip + np.ceil(trip['d'].iloc[-1]) + 1

plt.figure(figsize=[6, 2.5])
plt.plot(pair['d'], pair['eci'], 'r+')
if is_trip:
    plt.plot(trip['d'] + orig_for_trip, trip['eci'], 'b+')
if is_quad:
    plt.plot(quad['d'] + orig_for_quad, quad['eci'], 'g+')

ax = plt.gca()
locs = ax.xaxis.get_ticklocs()
locs = locs.tolist()

last_n = 0
for i in range(len(locs)):
    if is_trip and orig_for_trip > last_n and orig_for_trip < locs[i]:
        locs.insert(i, orig_for_trip)
        i += 1
    if is_quad and orig_for_quad > last_n and orig_for_quad < locs[i]:
        locs.insert(i, orig_for_quad)
        i += 1
    last_n = locs[i]

locs_labels = np.array(locs.copy())

if is_trip:
    locs_labels[locs > orig_for_trip] -= orig_for_trip
if is_quad:
    locs_labels[locs > orig_for_quad] -= (orig_for_quad - orig_for_trip)

locs_labels = locs_labels.astype('str')
locs_labels[0] = '\npair'
if is_trip:
    locs_labels[np.argmin(np.abs(locs - orig_for_trip))] = '\ntrip'
if is_quad:
    locs_labels[np.argmin(np.abs(locs - orig_for_quad))] = '\nquad'

ax.xaxis.set_ticks(locs)
ax.xaxis.set_ticklabels(locs_labels)

plt.axhline(0, 0, 1, linestyle=':', color='k')
plt.xlabel(r'Diameter ($\AA$)')
plt.ylabel('ECI (meV)')
plt.tight_layout()

if args.s:
    plt.savefig(args.s)
else:
    plt.show()
