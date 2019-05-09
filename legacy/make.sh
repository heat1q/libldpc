#!/bin/bash

name="sim_legacy"
defines=""
name="$1"
defines="$2 $3 $4 $5 $6 $7 $8 $9"
g++ -std=c++11 sim_legacy.cpp ldpc/ldpc.cpp ldpc/decoder.cpp sim/sim.cpp sim/startsim.cpp -o "$1" -O3 -w -pthread $defines
