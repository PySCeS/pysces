set(f2py_module_name "nleq2")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/nleq2.pyf"
                     "${CMAKE_CURRENT_SOURCE_DIR}/nleq2.f"
                     "${CMAKE_CURRENT_SOURCE_DIR}/linalg_nleq2.f"
                     "${CMAKE_CURRENT_SOURCE_DIR}/zibmon.f"
                     "${CMAKE_CURRENT_SOURCE_DIR}/zibsec.f"
                     "${CMAKE_CURRENT_SOURCE_DIR}/zibconst.f"
                     "${CMAKE_CURRENT_SOURCE_DIR}/wnorm.f")

set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
set(generated_libs_dir "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}/.libs")

# work around path expansion bug in Windows
# see https://groups.google.com/g/scikit-build/c/L7SmriBFeFk
# set(f2py_build_dir f2py_build)

add_custom_target(${f2py_module_name} ALL
  DEPENDS ${generated_module_file}
  )

add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND ${F2PY_EXECUTABLE}
    -c ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

install(FILES ${generated_module_file} DESTINATION "pysces/nleq2")
if (WIN32)
  install(DIRECTORY ${generated_libs_dir} DESTINATION pysces
        FILES_MATCHING PATTERN "*.dll")
endif()
