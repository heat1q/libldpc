#!/bin/bash

name="sim_cpu"
defines=""
name="$1"
defines="$2 $3 $4 $5 $6 $7 $8 $9"
g++ -std=c++11 sim_cpu.cpp ldpc/ldpc.cpp ldpc/decoder.cpp sim/ldpcsim.cpp sim/start.cpp -o $name -O3 -w -pthread $defines
