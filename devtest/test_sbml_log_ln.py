"""
See [https://github.com/PySCeS/pysces/issues/44] for details
"""
# input
model_file = 'BIOMD0000000012_url.xml'

# environment
import os
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import pysces

model_dir = os.path.join(cDir, 'models')
out_dir = os.path.join(cDir, 'out', model_file)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
pysces.output_dir = out_dir

# code
pysces.interface.convertSBML2PSC(model_file, sbmldir=model_dir, pscdir=out_dir)

mod = pysces.model(model_file+'.psc', dir=out_dir)
mod.doSim(1000, 1000)
mod.SimPlot(['PX','PY','PZ'])
pysces.plt.setRange('x', min=400, max=1000)
pysces.plt.export(model_file, directory=out_dir, type='png')

print('\nOutput path: \"{}\"\n\nDone.\n'.format(out_dir))

