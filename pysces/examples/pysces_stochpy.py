"""
Release 0.8.x of PySCeS includes a prototype interface to StochPy (stompy.sf.net)
developed by Timo Maarleveld at the VU University Amsterdam. Incorporating PySCeS technology,
StochPy adds a growing number of stochastic simulation and analysis capabilities to PySCeS. 
To activate the interface simply install StochPy, load a model suited for discrete simulation 
and the following methods become available on the PySCeS model object:

mod.doStochSim()
mod.StochSimPlot()
mod.doStochSimPlot(

The following example highlights some of the StochPy functionality available in PySCeS and 
uses the input file Burstmodel.psc (also available in the PySCeS example files). All the 
graphs generated will be both displayed and saved as PNG format files. 


NO WARRANTY OF FITNESS OF ANY KIND FOR ANYTHING ... SERIOUSLY!
(C) Brett G. Olivier Amsterdam 2011. All rights reserved
"""

# First we need to import PySCeS 0.8.0 or newer and store the current directory in CDIR
import os, time 
CDIR = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import numpy, pysces

# Next we instantiate a PySCeS model object (load the input file from the current directory)
# and set the simulation endtime/endsteps. Note how in this input file we specify that this 
# model should be simulated stochastically and that the output should be displayed in amounts. 
mod = pysces.model('Burstmodel.psc', CDIR)
MODEL_NAME = 'Burstmodel'
ENDTIME = 100
ENDSTEPS = 10000

# First let's assume this model represents a determinisitic system and do a continuous 
# simulation to endtime with 2*endtime time points with LSODA or CVODE. 
mod.doSimPlot(ENDTIME, 2*ENDTIME)
# Save the results to a file using the PySCeS plotting library (pysces.plt.*)
pysces.plt.export('%s_par_set_1_deterministic' % MODEL_NAME, CDIR)

# Typical continuous result. Now let's assume this is a stochastic system and do a discrete 
# simulation to the same end time using the Direct method 
mod.doStochSim(end=ENDTIME, mode='time', method='Direct')
# we can plot the results with this function and save the results
mod.StochSimPlot()
pysces.plt.export('%s_par_set_1_stochastic_endtime' % MODEL_NAME, CDIR)

# or we can do everything at once, lets try simulate to a number of steps and plot the results
mod.doStochSimPlot(end=ENDSTEPS, mode='steps', method='Direct')
pysces.plt.export('%s_par_set_1_stochastic_maxsteps' % MODEL_NAME, CDIR)

input('\npress <enter> to continue\n')

# a cool feature of PySCeS is a totally encapsulated model object which means
# we can create an independent object from the same file, let's instantiate mod2
mod2 = pysces.model('Burstmodel.psc', CDIR)

# This time however we want to change some parameters ... in PySCeS this is
# done in 'real' time

# Set parameter set 2
mod2.koff = 50
mod2.kon = 50
mod2.ksyn = 80
mod2.kdeg = 2.5

# again lets to a continous simulation to endtime with 2*endtime time points
# and save the results to a file
mod2.doSimPlot(ENDTIME, 2*ENDTIME)
pysces.plt.export('%s_par_set_2_deterministic' % MODEL_NAME, CDIR)

# However as a continous model we also have access to PySCeS other functionality, 
# like steady-state solution and structural properties. Note, this might not makes sense 
# in the context the model was built for but is used here for illustration. 
mod2.doStateShow()

# structural properties
mod2.showFluxRelationships()
mod2.showConserved()
mod2.showODEr()

# Metabolic Control Analysis
mod2.doMca()
mod2.showElas()
mod2.showCC()

input('\npress <enter> to continue\n')

# now lets switch back to a discrete simulation using the Direct method
mod2.doStochSimPlot(end=ENDTIME, mode='time', method='Direct')
pysces.plt.export('%s_par_set_2_stochastic_endtime' % MODEL_NAME, CDIR)

# of course we can also plot other interesting things, however, 
# now we use the StochSimPlot fucntion as we don't want to rerun the simulation
# simply plot different things like: waiting times
mod2.StochSimPlot(plot='waiting_times',format='points')
pysces.plt.export('%s_par_set_2_stochastic_waiting_times' % MODEL_NAME, CDIR)

# reaction propensities
mod2.StochSimPlot(plot='propensities',format='*')
pysces.plt.export('%s_par_set_2_stochastic_propensities' % MODEL_NAME, CDIR)

# Anything we like (here we plot a concentration and propensity) and format the
# resulting graph using the PySCeS plotting library: 
mod2.StochSimPlot(plot=['mRNA', 'R3'], format='.')
pysces.plt.setGraphTitle('My new plot')
pysces.plt.setRange('x',40,60)
pysces.plt.export('%s_par_set_2_stochastic_selected' % MODEL_NAME, CDIR)

input('\npress <enter> to exit')
