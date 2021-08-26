#!/usr/bin/env python3
"""
Test steady state (HYBRD and NLEQ2), MCA and parameter scans
"""

# input
model_file = 'rohwer_sucrose2'

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
# steady state
mod = pysces.model(model_file+'.psc', dir=model_dir)
mod.doState()
ss_data = np.copy(mod.data_sstate.getAllStateData())

# mca
mod.doMca()
elas_data = np.copy(mod.elas_var)
cc_data = np.copy(mod.cc_all)

# scan
mod.scan_in = 'Vmax9'
mod.scan_out = ['Suc_ss', 'Fru_ss', 'Glc_ss', 'J_R9', 'J_R11', 'ecR9_Suc',
                'ccJR9_R9', 'ccJR11_R9']
scan_range = np.linspace(0.1, 1.0, 91)
mod.Scan1(scan_range)
scan_data = np.copy(mod.scan_res)

# plots
mod.Scan1Plot(['Suc_ss', 'Fru_ss', 'Glc_ss'])
pysces.plt.export(model_file+'_species', directory=out_dir, outtype='png')
mod.Scan1Plot(['J_R9', 'J_R11'])
pysces.plt.export(model_file+'_fluxes', directory=out_dir, outtype='png')
mod.Scan1Plot(['ecR9_Suc', 'ccJR9_R9', 'ccJR11_R9'])
pysces.plt.export(model_file+'_mca', directory=out_dir, outtype='png')

ref_ss_data = pickle.load(open(data_dir+'/ssdata.pkl', 'rb'))
ref_elas_data, ref_cc_data = pickle.load(open(data_dir+'/mcadata.pkl', 'rb'))
ref_scan_data = pickle.load(open(data_dir+'/scandata.pkl', 'rb'))

assert np.allclose(ss_data, ref_ss_data), \
    "Steady-state data doesn't match reference!"
assert np.allclose(elas_data, ref_elas_data), \
    "Elasticity data doesn't match reference!"
assert np.allclose(cc_data, ref_cc_data), \
    "Control coefficient data doesn't match reference!"
assert np.allclose(scan_data, ref_scan_data), \
    "Scan data doesn't match reference!"

print('\nOutput path: \"{}\"\n\nDone.\n'.format(out_dir))
