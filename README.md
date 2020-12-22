# LDPC Code Simulation

![Build Status](https://github.com/heat1q/libldpc/workflows/Build/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

### Requirements
* CMake Version 3.10
* C++17 

### Getting Started
Building the binary:

`$ cmake . `
 
`$ make ldpc_sim` produces an executeable containing the simulator, which is used by e.g. the commandline. See **Running the Simulator**.

`$ make ldpc` produces a shared library containing the simulator for external usage. See **Python Wrapper**.

### Running the Simulator
After successful build the simulator can be executed. Note the usage:
```
$./ldpc_sim --help
Usage: ldpc [options] codefile output-file snr-range 

Positional arguments:
codefile            	LDPC codefile containing all non-zero entries, compressed sparse row (CSR) format.
output-file         	Results output file.
snr-range           	{MIN} {MAX} {STEP}

Optional arguments:
-h --help           	shows help message and exits
-v --version        	prints version information and exits
-G --gen-matrix     	Generator matrix file, compressed sparse row (CSR) format.
-i --num-iterations 	Number of iterations for decoding. (Default: 50)
-s --seed           	RNG seed. (Default: 0)
-t --num-threads    	Number of frames to be decoded in parallel. (Default: 1)
--channel           	Specifies channel: "AWGN", "BSC", "BEC" (Default: AWGN)
--decoding          	Specifies decoding algorithm: "BP", "BP_MS" (Default: BP)
--max-frames        	Limit number of decoded frames.
--frame-error-count 	Maximum frame errors for given simulation point.
--no-early-term     	Disable early termination for decoding.
```


### Python Wrapper
The simulator may be used as Python Module in a threaded application.
```
from pyLDPC import ldpc

# initialize with codefile & simfile
code = ldpc.LDPC("doc/exmp_code.txt")

# simulate AWGN with 100 iters and MinSum decoding
code.simulate(snr=[0, 6, 0.2], iterations=100, decoding="BP_MS")

# stop the simulation
code.stop_simulation()

# retrieve the results
print(code.get_results())
```
