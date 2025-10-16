project(
  'pysces',
  ['c', 'fortran'],
  version : run_command(
    ['pysces/version.py'],
    check: true
  ).stdout().strip(),
  meson_version: '>= 1.2.3'
)

fc = meson.get_compiler('fortran')

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(
  py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
fortranobject_c = incdir_f2py / 'fortranobject.c'

inc_np = include_directories(incdir_numpy, incdir_f2py)

numpy_nodepr_api = ['-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION']

fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  c_args: numpy_nodepr_api,
  dependencies: py_dep,
  include_directories: inc_np,
  gnu_symbol_visibility: 'hidden',
)

fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: inc_np,
)

fortran_modules = run_command(py,
  ['config.py'],
  check : true
).stdout()

subdir('pysces')
