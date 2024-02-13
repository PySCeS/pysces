nleq2_wrappers = custom_target(
  'nleq2module.c',
  input : ['nleq2.pyf'],
  output : ['nleq2module.c', 'nleq2-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@']
)

py.extension_module(
  'nleq2',
  [
    'nleq2.f',
    'linalg_nleq2.f',
    'zibmon.f',
    'zibsec.f',
    'zibconst.f',
    'wnorm.f',
    nleq2_wrappers,
    fortranobject_c
  ],
  include_directories: inc_np,
  dependencies : [
    py_dep,
    quadmath_dep,
  ],
  install : true,
  subdir: 'pysces/nleq2'
)

py.install_sources(
  [
    '__init__.py'
  ],
  subdir: 'pysces/nleq2'
)

install_data(
  ['nleq2_readme.txt', 'readme.txt'],
  install_dir: pysces_dir / 'nleq2'
)
