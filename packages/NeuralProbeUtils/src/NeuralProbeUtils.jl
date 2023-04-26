module NeuralProbeUtils

using DSP, Statistics#, Polyester
using DataChunks, TensorOps

export preprocessor, load_and_process, channel_order, find_bad_channels,
    interpolate_bad_channels!

include("./preprocessing.jl")
include("./bad_channels.jl")

abstract type ProbeData{T} end

# basic data dimentions
n_channel(p::ProbeData) = error("Not implemented")

# read access
memmap(p::ProbeData) = error("Not implemented")

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

    trial_groups = group(idx, out.sizes)

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

    #Threads.@threads
    Threads.@threads for k in 1:nchunks(out)
        for (j, slice) in enumerate(eachslice(getchunk(out,k), dims=3))
            # TODO we need padding....
            start = trial_groups[k][j] - npre
            stop = start + len - 1
            preprocess!(slice, raw[korder,start:stop]', procs[k])
        end
    end

    return out
end

function group(x::AbstractVector{T}, grps::AbstractVector{<:Integer}) where {T}
    out = Vector{Vector{T}}(undef, length(grps))
    j = 1
    for k in 1:length(out)
        out[k] = x[j:grps[k]]
        j = grps[k] + 1
    end
    return out
end

end
