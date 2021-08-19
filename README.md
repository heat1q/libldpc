# LDPC Code Simulation

![Build Status](https://github.com/heat1q/libldpc/workflows/Build/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

### Requirements
* CMake Version 3.10
* C++17 

### Getting Started
Building the binary:

```
$ cmake .
$ cmake --build . --target TARGET
```
where the target can be the following:
* `--target ldpcsim` produces an executeable containing the simulator, which is used by e.g. the commandline. See **Running the Simulator**.

* `--target ldpc` produces a shared library containing the simulator for external usage. See **Python Wrapper**.

### Running the Simulator
After successful build the simulator can be executed. Note the usage:
```
$ ./ldpcsim --help
Usage: ldpc [options] codefile output-file snr-range 

Positional arguments:
codefile                LDPC parity-check matrix file containing all non-zero entries.
output-file             Results output file.
snr-range               {MIN} {MAX} {STEP}

Optional arguments:
-h --help               shows help message and exits
-v --version            prints version information and exits
-G --gen-matrix         Generator matrix file.
-i --num-iterations     Number of iterations for decoding. (Default: 50)
-s --seed               RNG seed. (Default: 0)
-t --num-threads        Number of frames to be decoded in parallel. (Default: 1)
--channel               Specifies channel: "AWGN", "BSC", "BEC" (Default: AWGN)
--decoding              Specifies decoding algorithm: "BP", "BP_MS" (Default: BP)
--max-frames            Limit number of decoded frames.
--frame-error-count     Maximum frame errors for given simulation point.
--no-early-term         Disable early termination for decoding.
```

### Code File Format
The code file lists the the non-zero entries of the matrix line-by-line, 
where the left number refers to the row index and the right to the column index. For example, the matrix 
```
0 0 1 0
0 1 0 0
0 0 0 1
1 1 1 0
```
is represented as
```
0 2
1 1
2 3
3 0
3 1
3 2
```
Optionally, a puncturing (or shortening) patter may be specified a the top of the file:

*Puncture 2 bits associated with column index 1 and 3*
```
puncture [2]: 1 3 
```
A sample `k=128`, `n=1024` LDPC code with parity-check matrix `h.txt` and generator matrix `g.txt` can be found in `tests/code/`.


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

### Contributing
If you found a bug or have any suggestions or improvements, feel free to start a discussion or submit a PR!