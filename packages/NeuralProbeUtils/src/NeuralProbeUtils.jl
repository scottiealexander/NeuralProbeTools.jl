module NeuralProbeUtils

using DSP, Statistics, Mmap#, Polyester
using DataChunks, TensorOps

import Base

const RealMat = AbstractMatrix{<:Real}

export FlatBinaryFile, Neuropixel3A, DBCDeepArray, channel_order, vertical_pitch,
    reference_channels, channel_positions

export AbstractProcessor, Preprocessor, Resampler, Filterer, Interpolator, DepthAverager

export load_and_process, find_bad_channels, interpolate_bad_channels!, fileindex_to_depthindex

include("./probe_definitions.jl")
include("./flatbinary_file.jl")
include("./processors.jl")
include("./preprocessing.jl")
include("./bad_channels.jl")

function load_and_process(d::FlatBinaryFile{P,T}, idx::AbstractVector{<:Integer},
    npre::Integer, npost::Integer, proc::AbstractProcessor) where {P,T}

    # size of time dimention in read-in data
    len = npost + npre

    # channel permutation vector
    korder = channel_order(d)

    olen, nchan = output_size(proc, len, length(korder))

    out = DataChunk(output_type(T, proc), (olen, nchan, length(idx)), ChunkDim(3))

    trial_groups = DataChunk(idx, ChunkDim(1), nchunks=nchunks(out))

    raw = memmap(d)

    # make nchunks(out) copies of the processor so that mutable data is not
    # shared between tasks
    procs = map(x -> copy(proc), 1:nchunks(out))

    Threads.@threads for k in 1:nchunks(out)
        chnk = getchunk(out, k)
        onset = getchunk(trial_groups, k)
        for (j, slice) in enumerate(eachslice(chnk, dims=3))
            start = onset[j] - npre
            stop = start + len - 1

            # yep, views of mmap arrays work just fine...
            process!(slice, view(raw, korder, start:stop)', procs[k])
        end
    end

    return out
end

#=
function thread_error_wme()

    x = DataChunk(1:27, ChunkDim(1))

    # # this causes the error, where the name <ychnk> is repeated here and in the threaded
    # # loop (NOTE: putting local ychnk = ... with the loop also avoid the error)
    ychnk = 1:27
    y = DataChunk(ychnk, ChunkDim(1))

    # # this does not
    # y = DataChunk(1:27, ChunkDim(1))

    @assert(nchunks(x) == nchunks(y))
    @assert(x == y)

    Threads.@threads for k in 1:nchunks(x)
        xchnk = getchunk(x, k)
        # local ychnk = getchunk(y, k) # <- use of local avoids error...
        ychnk = getchunk(y, k)
        for j in eachindex(xchnk)
            @assert(xchnk[j] == ychnk[j], "$(xchnk[j]) != $(ychnk[j])")
        end
    end
    return nothing
end
=#
end
