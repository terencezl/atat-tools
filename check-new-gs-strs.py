#!/bin/env python

import os
import subprocess
import argparse
import numpy as np
import pandas as pd
import atat_module
from scipy import interpolate
from scipy.spatial import ConvexHull

parser = argparse.ArgumentParser(
    description="Prepare the directories for the new ground states to run VASP.")
parser.add_argument('--genstr', default='genstr.out', help="the genstr.out file")
parser.add_argument('-x', action='store_true', help="actually create the directories.")
args = parser.parse_args()

strcorrenergy = pd.read_csv('strcorrenergy.csv', index_col='str_id')
gs = atat_module.get_df('gs', 'gs.out')
interp_func = interpolate.interp1d(gs['c'], gs['E_fit'])
strcorrenergy['E_gs'] = interp_func(strcorrenergy['c'])
strcorrenergy['E_diff'] = strcorrenergy['E_per_site'] - strcorrenergy['E_gs']
str_below_hull = strcorrenergy[(strcorrenergy['E_diff'] < -0.0001)]

str_id_list = []
for c in sorted(list(set(str_below_hull['c']))):
    str_below_hull_c = str_below_hull[np.isclose(str_below_hull['c'], c)]
    str_id_list.append(str_below_hull_c['E_diff'].idxmin())

if len(str_id_list) >= 3:
    str_lowest_hanging = strcorrenergy.loc[str_id_list]
    str_gs = gs[['c']].copy()
    str_gs['E_diff'] = 0
    concat_df = pd.concat([str_lowest_hanging, str_gs])
    hull = ConvexHull(concat_df[['c', 'E_diff']].values)
    str_id_list = concat_df.iloc[hull.vertices].index.tolist()

for str_id in gs.index:
    if str_id in str_id_list:
        str_id_list.remove(str_id)

new_gs_df = str_below_hull.loc[str_id_list]
print(new_gs_df)

if args.x:
    counter = 0
    with open(args.genstr) as f:
        for line in f:
            if counter in str_id_list:
                line = next(f)
              #  dirname = 'gs' + str(counter)
                dirname = str(counter)
                if not os.path.isdir(dirname):
                    os.mkdir(dirname)
                    fout = open(os.path.join(dirname, 'str.out'), 'w')
                    while 'end' not in line:
                        fout.write(line)
                        line = next(f)
                    fout.close()
                    subprocess.getoutput('touch ' + os.path.join(dirname, 'wait'))
                else:
                    while 'end' not in line:
                        line = next(f)
            if 'end' in line:
                counter += 1
