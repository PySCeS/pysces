pitcon_wrappers = custom_target(
  'pitconmodule.c',
  input : ['pitcon.pyf'],
  output : ['pitconmodule.c', 'pitcon-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@']
)

py.extension_module(
  'pitcon',
  [
    'pcon61subd.f'
    'dpcon61.f'
    'dpcon61w.f'
    pitcon_wrappers,
    fortranobject_c
  ],
  include_directories: inc_np,
  dependencies : [
    py_dep,
    quadmath_dep,
  ],
  install : true,
  subdir: 'pysces/pitcon'
)

py.install_sources(
  [
    '__init__.py'
  ],
  subdir: 'pysces/pitcon'
)
