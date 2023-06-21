# NeuralProbeTools

A collection of Julia modules for accessing and analyzing high-density neural probe data, including a highly optimized and extensible processing pipeline framework.

## Organization
* `utils` - general purpose utilities (i.e. not specific to neural probe data) that included packages and analyses rely upon
* `packages` - core functionality for loading and processing neural probe data
    * `NeuralProbeUtils` - core routines applicable to all probe data (i.e. any probe + file format combination) including:
        * Probe definitions
        * Flat binary file interface
        * Processor / Preprocessing interface and implementations
    * `GLX` - Access datasets recorded by SpikeGLX
    * `Intan` - Access Intan datasets
    * `OEphys` - Access OpenEphys datasets
* `analysis` - User facing modules for performing analyses (see function level help within modules for more documentation)
    * `OEphysCSD` - calculate and display current source density (CSD) images for (e.g.) during-acquisition probe positioning


## Install notes

```bash
# repo URL
git clone https://github.com/scottiealexander/NeuralProbeTools.jl NeuralProbeTools

# launch julia (e.g.) 1.8.5 with (e.g.) 4 threads
julia-1.8.5 -t 4
```

Then in julia:

```julia
# where ever you installed the NeuralProbeTools repo
neural_tools_dir = "..."

# if the OEphysCSD package has *NOT* been compiled into the sysimage, these
# commands are needed (each time) after launching a new Julia session
cd(neural_tools_dir)

# add project code to julia's search path (needed every time)
include("./setpath.jl")

# ---------------------------------------------------------------------------- #
# actually download and install dependencies (ONLY NEEDED THE FIRST TIME)
import Pkg
Pkg.instantiate()
# ---------------------------------------------------------------------------- #

# pull in main CSD code (needed every time)
using OEphysCSD

# now we can run the code as many times as we need / like
rec_dir = "..." # full path to an OpenEphys data directory

# outputs:
#   h -> figure handle
#   erp -> time x channels x trials data array
#   bad -> vector of indices of "bad" channels
h, erp, bad = OEphysCSD.run(rec_dir)

# or type '?' as the first character into the repl to get a "help prompt", then type the function name `OEphysCSD.run` and [ENTER]

# to time code (prints runtime and allocation info) you can use:
@time h, erp, bad = OEphysCSD.run(rec_dir)

```

## Notes

### Bad channels / channel indices
With in a processing pipeline, channel indices range from `1:n_channel`, where channel 1 is farthest from the probe tip, and channel `n_channel` is closest
to the probe tip, as defined in `NeuralProbeUtils/src/probe_definitions.jl`.

So to convert from manufacturer / data file indices (assuming 1-based indexing everywhere) to "depth indices" you would:

```julia
using NeuralProbeUtils
probe = DBCDeepArray() # or whatever you're using

# for a DBCDeepArray, in the data file the channels are stored such that
# channel 1 is farthest from probe tip, 128 is second farthest, 65 is closest, etc.
file_idx = [1, 128, 64, 65]

depth_idx = fileindex_to_depthindex(probe, file_idx) # -> [1, 2, 127, 128]

```

### Custom sysimage

Launch julia, then:

```julia
import Pkg #(or type ']' as the first character into the repl to get the Pkg prompt and then omit the 'Pkg.' from what follows)
Pkg.add("PackageCompiler")

# then as above:
cd(neural_tools_dir)
include("./setpath.jl") # make sure PackageCompiler can locate the OEphysCSD package for the next line

#then run PackageCompiler
PackageCompiler.create_sysimage(["OEphysCSD"]; # name of "main" Package to compile
    sysimage_path="../test_sysimage.so", # output path for the custom sysimage (perhaps change extension to dll for Windows)
    precompile_execution_file="./precompile.jl" # script from which to build a list of functions / Packages to compile into the sysimage
)
```

then exit julia and relaunch:

```bash
julia-1.8.5 -t 4 -J../test_sysimage.so # or whatever path and extension
```

## To do
* Padding for in-place filtering during preprocessing
* ~~Allow for custom Preprocessor pipelines?~~
* ~~Add utils for spike glx data sets~~
* ~~Separate acquisition / file format code from probe code?~~
* ~~Automatic detection of "bad" channels?~~
