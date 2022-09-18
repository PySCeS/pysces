set(f2py_module_name "pitcon")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/pitcon.pyf"
                     "${CMAKE_CURRENT_SOURCE_DIR}/pcon61subd.f"
                     "${CMAKE_CURRENT_SOURCE_DIR}/dpcon61.f"
                     "${CMAKE_CURRENT_SOURCE_DIR}/dpcon61w.f")

set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
set(generated_libs_dir "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}/.libs")

# work around path expansion bug in Windows
# see https://groups.google.com/g/scikit-build/c/L7SmriBFeFk
set(f2py_build_dir f2py_build)

add_custom_target(${f2py_module_name} ALL
  DEPENDS ${generated_module_file}
  )

add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    -c ${fortran_src_file}
    --build-dir ${f2py_build_dir}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

install(FILES ${generated_module_file} DESTINATION "pysces/pitcon")
if (WIN32)
  install(DIRECTORY ${generated_libs_dir} DESTINATION pysces
        FILES_MATCHING PATTERN "*.dll")
endif()
