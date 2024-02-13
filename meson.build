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
quadmath_dep = fc.find_library('quadmath', required: false)

fortran_modules = run_command(py,
  ['config.py'],
  check : true
).stdout()

subdir('pysces')
