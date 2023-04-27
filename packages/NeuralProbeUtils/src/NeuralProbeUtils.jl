module NeuralProbeUtils

using DSP, Statistics, Mmap#, Polyester
using DataChunks, TensorOps

export preprocessor, load_and_process, channel_order, find_bad_channels,
    interpolate_bad_channels!

include("./preprocessing.jl")
include("./bad_channels.jl")

abstract type ProbeData{T} end

# basic data dimentions
n_channel(p::ProbeData) = error("n_channel() not implemented")

# read access
memmap(p::ProbeData) = error("memmap() not implemented")

# channel order / reordering
channel_order(p::ProbeData) = 1:n_channel(p)

function load_and_process(d::ProbeData{T}, idx::AbstractVector{<:Integer},
    npre::Integer, npost::Integer, proc::Preprocessor=preprocessor()) where {T}

    # size of time dimention in read-in data
    len = npost + npre

    # channel permutation vector
    korder = channel_order(d)

    olen = output_length(proc, len)
    L = output_type(T, proc)

    out = DataChunk(L, (olen, length(korder), length(idx)), ChunkDim(3))

    trial_groups = DataChunk(idx, ChunkDim(1), nchunks=nchunks(out))

    raw = memmap(d)

    # FIRFilter within the Resampler (within <proc>) is mutable, so each
    # thread needs their own, others can be references to a single instance
    procs = map(1:nchunks(out)) do _
        Preprocessor(
            resampler=deepcopy(proc.resampler),
            filter=proc.filter,
            interpolator=proc.interpolator
            )
    end

    Threads.@threads for k in 1:nchunks(out)
        chnk = getchunk(out, k)
        onset = getchunk(trial_groups, k)
        for (j, slice) in enumerate(eachslice(chnk, dims=3))
            start = onset[j] - npre
            stop = start + len - 1
            # yep, views of mmap arrays work just fine...
            preprocess!(slice, view(raw, korder, start:stop)', procs[k])
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
