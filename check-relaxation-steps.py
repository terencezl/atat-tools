#!/bin/env python
import sys
import pandas as pd
from glob import glob
import pymatgen as mg
import subprocess
import argparse

parser = argparse.ArgumentParser(
    description="Check what structures have relaxation steps larger or equal to a certain integer.")
parser.add_argument('-n', type=int, help="number of steps, default to the maximum.")
args = parser.parse_args()

result_list = []

for i in glob('*/OSZICAR.relax'):
    result_list.append([int(i.strip('/OSZICAR.relax')), len(mg.io.vasp.Oszicar(i).ionic_steps)])

df = pd.DataFrame(result_list, columns=['str_id', 'n_ionic_steps']).set_index('str_id').sort_index()
print(df.describe())
if args.n:
    print('Structures with the specified number, ' + args.n + ', of relaxation steps:')
    indices = df[df['n_ionic_steps'] >= int(args.n)].index.values
else:
    print('Structures with the most relaxation steps:')
    indices = df[df['n_ionic_steps'] == df['n_ionic_steps'].max()].index.values

print(indices)

for i in indices:
    print(i)
    subprocess.call("grad2 " + str(i) + '/OUTCAR.relax.gz | tail -3', shell=True)
