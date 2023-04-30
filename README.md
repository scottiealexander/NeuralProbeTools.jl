# To do
* Padding for in-place filtering during preprocessing
* Allow for custom Preprocessor pipelines?
* Add utils for spike glx data sets
* ~~Separate acquisition / file format code from probe code?~~
* ~~Automatic detection of "bad" channels?~~

```julia
# ============================================================================ #
function spikeglx_data(basedir::AbstractString, ap::Bool=false)
    # ...
    return FlatBinaryFile{Neuropixel3A,Int16}(Neuropixel3A(), filepath, nchan, nsample, fs)
end
# ============================================================================ #
```
