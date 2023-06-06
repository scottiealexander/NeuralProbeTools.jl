## To do
* Padding for in-place filtering during preprocessing
* ~~Allow for custom Preprocessor pipelines?~~
* ~~Add utils for spike glx data sets~~
* ~~Separate acquisition / file format code from probe code?~~
* ~~Automatic detection of "bad" channels?~~

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

include("./setpath.jl")
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
