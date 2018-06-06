#!/bin/env python
import pandas as pd
from subprocess import getoutput

natoms_per_unitcell = int(getoutput('cellcvrt -pn < lat.in'))
nsites_per_unitcell = int(getoutput('grep , lat.in | wc -l'))
divisor = natoms_per_unitcell / nsites_per_unitcell

strcorr = pd.read_table('strcorr.out', header=None, sep='\s+')
strenergy = pd.read_table('strenergy.out', header=None, sep='\s+')
natoms = pd.read_table('natoms.out', header=None, sep='\s+')
strcorrenergy = pd.concat([(strcorr[1] + 1) / 2, strenergy[0], natoms[0]], 1, ignore_index=True)
strcorrenergy.columns = ['c', 'E', 'natoms']
strcorrenergy['nsites'] = strcorrenergy['natoms']/divisor
strcorrenergy['E_per_site'] = strcorrenergy['E'] / nsites_per_unitcell
strcorrenergy.index.name = 'str_id'
strcorrenergy.to_csv('strcorrenergy.csv')
