#!/usr/bin/env python3
"""
See [https://github.com/PySCeS/pysces/issues/44] for details
"""

# input
model_file = 'BIOMD0000000012_url.xml'

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
# explicitly call LSODA as CVODE gets invoked as soon as there are assignment rules
mod.mode_integrator = 'LSODA'
mod.doSim(1000, 1000)
lsoda_data = np.copy(mod.data_sim.getAllSimData())
mod.SimPlot(['PX','PY','PZ'])
pysces.plt.setRange('x', min=400, max=1000)
pysces.plt.export(model_file, directory=out_dir, outtype='png')

mod.reLoad()
mod.mode_integrator='CVODE'
mod.__settings__['cvode_track_assignment_rules'] = False
mod.doSim(1000, 1000)
cvode_data = np.copy(mod.data_sim.getAllSimData())

ref_sim_data = pickle.load(open(data_dir+'/simdata.pkl', 'rb'))
assert np.allclose(lsoda_data, ref_sim_data), \
    "LSODA data doesn't match reference!"
assert np.allclose(cvode_data, ref_sim_data), \
    "CVODE data doesn't match reference!"

print('\nOutput path: \"{}\"\n\nDone.\n'.format(out_dir))

