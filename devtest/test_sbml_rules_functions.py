#!/usr/bin/env python3
"""
Test SBML assignment rules and functions
"""

# input
model_file = 'BIOMD0000000860.xml'

# environment
import os
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import pysces
import numpy as np
import pickle

model_dir = os.path.join(cDir, 'models')
out_dir = os.path.join(cDir, 'out', model_file)
data_dir = os.path.join(cDir, 'ref_data', model_file)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
pysces.output_dir = out_dir

# code
pysces.interface.convertSBML2PSC(model_file, sbmldir=model_dir, pscdir=out_dir)

mod = pysces.model(model_file+'.psc', dir=out_dir)
mod.doSim(4*3600, 4*3600+1)
mod.SimPlot()
pysces.plt.export(model_file, directory=out_dir, outtype='png')

ref_sim_data = pickle.load(open(data_dir+'/simdata.pkl', 'rb'))
assert np.allclose(mod.data_sim.getAllSimData(), ref_sim_data), \
    "Data doesn't match reference!"

print('\nOutput path: \"{}\"\n\nDone.\n'.format(out_dir))
