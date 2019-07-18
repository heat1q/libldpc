# Cuda LDPC Code Simulation
### Structure
* `/legacy` Simulator for decoding on CPU (decode in layers optional)
* `/src` Simulator for decoding on GPU (decode in layers optional)

### Compiling
*Compiling steps for `/legacy` are identical*

`$sh src/make.sh NAME FLAGS`

`NAME` is the name of the gernerated executable.

`FLAGS` are compile time constants. Note that the prefix `-D` has to be added to flags, in order for the compiler to recognize them. E.g. compiling with the `GPU_ID` and `LOG_TP` flag set to `1`, we write

`$sh src/make.sh sim "-DGPU_ID=1 -DLOG_TP=1"`

All possible flags are listed in `flags.txt` with their respective default value and short desciption.

### Running
After successful build the simulator can be executed by 

`$./src/NAME -code CODEFILE -sim SIMFILE -map MAPPINGFILE -layer LAYERFILE -threads NUMTHREADS`

The arguments `-layer` and `-threads` are optional.

`LAYERFILE` defines the subsets of check nodes that form a layer. A template can be found in `codes/exmp_layer.txt`
`NUMTHREADS` defines the number of frames processed in parallel (Not available for CPU decoding!)
