cmake_minimum_required(VERSION 3.10)

# set the project name
project(LDPC-Simulator)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -Wall -pthread")

set(SIM_FLAGS "" CACHE STRING "Compile flags (see flags.txt)")

# add the executable
add_executable(ldpc_sim src/sim_cpu.cpp src/core/functions.cpp src/core/ldpc.cpp src/decoding/decoder.cpp src/sim/ldpcsim.cpp src/sim/start.cpp)
target_compile_definitions(ldpc_sim PRIVATE LOG_FRAME_TIME=1 ${SIM_FLAGS})

# shared library
add_library(ldpc SHARED src/shared.cpp src/core/functions.cpp src/core/ldpc.cpp src/decoding/decoder.cpp src/sim/ldpcsim.cpp src/sim/start.cpp)
target_compile_definitions(ldpc PRIVATE LOG_FRAME_TIME=1 LIB_SHARED=1 ${SIM_FLAGS})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
