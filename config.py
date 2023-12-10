import os
import sys
import configparser

##### USER CONFIGURATION SECTION #####
# Set any of the below to False to disable the build and install of
# the corresponding module.
nleq2 = True
pitcon = True
##### END USER CONFIGURATION SECTION #####

if 'PYSCES_SKIP' in os.environ:
    if 'nleq2' in os.environ['PYSCES_SKIP']:
        nleq2 = False
    if 'pitcon' in os.environ['PYSCES_SKIP']:
        pitcon = False

local_path = os.path.dirname(os.path.abspath(__file__))
os.chdir(local_path)

# Default configurations for the pyscfg.ini files
config = {
    "gnuplot_dir": None,
    "silentstart": False,
    "change_dir_on_start": False,
    "custom_datatype": None,
}

def writeConfig(local_path, config={}):
    output = ''
    cp = configparser.ConfigParser()
    # PySCeS internal setup
    cp.add_section('Pysces')
    for key in config:
        print(repr(key) + ' :: ' + str(config[key]), file=sys.stderr)
        cp.set('Pysces', key, str(config[key]))
    # add configuration data
    cp.add_section('PyscesConfig')
    cp.set('PyscesConfig', 'matplotlib', 'True')
    # OSX patch thanks to AF
    if os.sys.platform == 'darwin':
        cp.set('PyscesConfig', 'matplotlib_backend', 'MacOSX')
    else:
        cp.set('PyscesConfig', 'matplotlib_backend', 'TkAgg')
    cp.set('PyscesConfig', 'gnuplot', 'False')
    # Built in modules
    cp.add_section('PyscesModules')
    if not pitcon:
        cp.set('PyscesModules', 'pitcon', 'False')
    else:
        cp.set('PyscesModules', 'pitcon', 'True')
        output += 'pitcon '
    # PySCeS external module setup
    cp.add_section('ExternalModules')
    if not nleq2:
        cp.set('ExternalModules', 'nleq2', 'False')
    else:
        cp.set('ExternalModules', 'nleq2', 'True')
        output += 'nleq2 '
    with open(os.path.join(local_path, 'pysces', 'pyscfg.ini'), 'w') as cfgfile:
        cp.write(cfgfile)
    return output

output = writeConfig(local_path, config)
print('Default configuration file installed', file=sys.stderr)
print(output)
