pysces_dir = py.get_install_dir() / 'pysces'

py.install_sources([
    'PyscesConfig.py',
    'PyscesContrib.py',
    'PyscesContribUser.py',
    'PyscesInterfaces.py',
    'PyscesJWSParse.py',
    'PyscesLink.py',
    'PyscesMiniModel.py',
    'PyscesModel.py',
    'PyscesModelMap.py',
    'PyscesParScan.py',
    'PyscesParse.py',
    'PyscesPlot.py',
    'PyscesPlot2.py',
    'PyscesRandom.py',
    'PyscesSBML.py',
    'PyscesSED.py',
    'PyscesScan.py',
    'PyscesStoich.py',
    'PyscesTest.py',
    'PyscesUtils.py',
    'PyscesWeb.py',
    'RateChar.py',
    '__init__.py',
    '_multicorescan.py',
    'version.py',
  ],
  subdir: 'pysces'
)

install_data(
  ['pyscfg.ini', 'version.txt'],
  install_dir: pysces_dir
)

pure_subdirs = [
  'contrib',
  'core2',
  'examples',
  'lib',
  'kraken',
  'pscmodels',
  'sandbox',
  'tests',
]

foreach subdir: pure_subdirs
  install_subdir(
    subdir,
    install_dir: pysces_dir,
    exclude_directories: ['sbw', 'visualise', '__pycache__']
  )
endforeach

# documentation
install_data(
  'docs/userguide.pdf',
  install_dir: pysces_dir / 'docs'
)

# nleq2
if fortran_modules.contains('nleq2')
  subdir('nleq2')
endif

# pitcon
if fortran_modules.contains('pitcon')
  subdir('pitcon')
endif
