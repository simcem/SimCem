cmake_minimum_required(VERSION 2.8.12)
project(simcem)
add_subdirectory(src/pybind11)
pybind11_add_module(core pysrc/simcem/core.cpp)
