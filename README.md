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
After successful build the simulator can be executed by  
`$./ldpc_sim -code [CODEFILE] -sim [SIMFILE] [-threads [NUMTHREADS] -seed [SEED]]`  

The arguments `-threads` and `-seed` are optional.\
`CODEFILE` Example codefile can be found in `doc/example_codefile.txt`\
`SIMFILE` Example simfile can be found in `doc/example_simfile.txt` \
`SEED` Random Number Generator (RNG) starting seed\
`NUMTHREADS` Number of frames processed in parallel\

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


