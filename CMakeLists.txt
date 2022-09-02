cmake_minimum_required(VERSION 3.17.0)

project(pysces)

find_package(PythonExtensions REQUIRED)
find_package(NumPy)
find_package(F2PY REQUIRED)

message(STATUS "F2PY included from: ${F2PY_INCLUDE_DIR}")

add_subdirectory("pysces/nleq2")
add_subdirectory("pysces/pitcon")
