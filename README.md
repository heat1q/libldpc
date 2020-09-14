# LDPC Code Simulation
### Structure
* `/gpu` Simulator for decoding on GPU (Not maintained currently)
* `/src` Simulator for decoding on CPU

### Compiling

`$ cmake . [-DSIM_FLAGS="[FLAGS]"]`
 
`FLAGS` are compile time constants. Multiple flags are to be separated by space and explicitly set to a value if used. E.g. `$ cmake . -DSIM_FLAGS="CN_APPROX_MINSUM=1 LOG_TP=1"`

All possible flags are listed in `flags.txt` with their respective default value and short desciption.

`$ make ldpc_sim` produces an executeable containing the simulator, which is used by e.g. the commandline. See a)

`$ make ldpc` produces a shared library containing the simulator for external usage. See b)

### a) Running the Executable
After successful build the simulator can be executed. Note the usage:
```
$./ldpc_sim --help
Usage: ldpc_sim [options] codefile output-file snr-range 

Positional arguments:
codefile            	LDPC codefile containing all non-zero entries, column major ordering.
output-file         	Results output file.
snr-range           	{MIN} {MAX} {STEP}

Optional arguments:
-h --help           	shows help message and exits
-v --version        	prints version information and exits
-i --num-iterations 	Number of iterations for decoding. (Default: 50)
-s --seed           	RNG seed. (Default: 0)
-t --num-threads    	Number of frames to be decoded in parallel. (Default: 1)
--channel           	Specifies channel {AWGN, BSC, BEC}
--decoding          	Specifies decoding algorithm {BP,HD}
--max-frames        	Limit number of decoded frames.
--frame-error-count 	Maximum frame errors for given simulation point.
--no-early-term     	Disable early termination for decoding.
```


### b) Use with Python threads
The simulator may be used as Python Module in a threaded application.
```
from pyLDPC import ldpc

# initialize with codefile & simfile
code = ldpc.LDPC("doc/exmp_code.txt", "doc/exmp_sim.txt")

# simulate in 1 Thread, set RNG seed to 0
code.simulate(1, 0)

# stop the simulation
code.stop_simulation()

# retrieve the results
print(code.get_results())
```


