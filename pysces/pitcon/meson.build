#pitcon_wrappers = custom_target(
#  'pitconmodule.c',
#  input : ['pitcon.pyf'],
#  output : ['pitconmodule.c', 'pitcon-f2pywrappers.f'],
#  command : [py, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@']
#)

pitcon_lib = static_library('pitcon_lib',
  'pcon61subd.f',
  'dpcon61.f',
  'dpcon61w.f',
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
)

py.extension_module(
  'pitcon',
  'pitconmodule.c',
  link_with: [pitcon_lib],
  dependencies: [py_dep, fortranobject_dep],
  install: true,
  link_language: 'fortran',
  include_directories: inc_np,
  subdir: 'pysces/pitcon'
)

py.install_sources(
  [
    '__init__.py'
  ],
  subdir: 'pysces/pitcon'
)
