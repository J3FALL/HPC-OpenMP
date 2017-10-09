cmake_minimum_required(VERSION 3.8)
project(HPC_OpenMP)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp)
add_executable(HPC_OpenMP ${SOURCE_FILES})