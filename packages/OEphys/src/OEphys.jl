module OEphys

using Mmap, DSP, JSON
using CSV, NeuralProbeUtils

export init_data, preprocessor, load_and_process

# ============================================================================ #
function openephys_data(basedir::AbstractString, d_codes::Dict{Int, String}, p::NeuralProbeUtils.AbstractProbe=DBCDeepArray())
    meta = metadata(basedir)
    filepath, nchan, nsample, fs, cont_idx = continuous_info(basedir, meta, 1)

    onset, dur, code, str = ttl_events_coded(basedir, d_codes, meta)
    #onset, dur, code = ttl_events(basedir, meta)
    onset .-= (cont_idx[1] - 1)

    return FlatBinaryFile{typeof(p),Int16}(p, filepath, nchan, nsample, fs), onset, dur, code, str
end
# ============================================================================ #
metadata(basedir::AbstractString) = JSON.parsefile(joinpath(basedir, "structure.oebin"))
# ---------------------------------------------------------------------------- #
function ttl_events(basedir::AbstractString, meta::Dict{String,Any}=metadata(basedir))

    ttlinfo = find_channel_by_name(meta["events"], "Rhythm FPGA TTL Input")
    ttldir = joinpath(basedir, "events", ttlinfo["folder_name"])

    sample_rate = ttlinfo["sample_rate"]::Float64

    states = read_npy_1D(Int16, joinpath(ttldir, "states.npy"))
    samples = read_npy_1D(Int64, joinpath(ttldir, "sample_numbers.npy"))
    # evt = read_npy_1D(Float64, joinpath(ttldir, "timestamps.npy"))

    nonset = sum(>(0), states)
    onset = Vector{Int}(undef, nonset)
    dur = Vector{Int}(undef, nonset)
    code = Vector{Int16}(undef, nonset)

    inc = 1
    for k in eachindex(states)
        if states[k] > 0
            kn = findnext(isequal(-states[k]), states, k+1)
            if kn != nothing
                onset[inc] = samples[k]
                dur[inc] = round(Int, ((samples[kn] - samples[k]) / sample_rate) * 1000)
                code[inc] = states[k]
                inc += 1
            end
        end
    end

    # sample onset (indices), pulse duration (in milliseconds), event code
    return onset, dur, code
end
#----------------------------------------------------------------------------- #
function ttl_events_coded(basedir::AbstractString, d::Dict{Int, String}, meta::Dict{String,Any}=metadata(basedir))

    ttlinfo = find_channel_by_name(meta["events"], "Rhythm FPGA TTL Input")
    ttldir = joinpath(basedir, "events", ttlinfo["folder_name"])

    sample_rate = ttlinfo["sample_rate"]::Float64

    states = read_npy_1D(Int16, joinpath(ttldir, "states.npy"))
    samples = read_npy_1D(Int64, joinpath(ttldir, "sample_numbers.npy"))
    # evt = read_npy_1D(Float64, joinpath(ttldir, "timestamps.npy"))


    nonset = sum(>(0), states)
    onset = Vector{Int}(undef, nonset)
    dur = Vector{Int}(undef, nonset)
    code = Vector{Int16}(undef, nonset)
    str = ""

    inc = 1
    for k in eachindex(states)
        if states[k] > 0
            kn = findnext(isequal(-states[k]), states, k+1)
            if kn != nothing
                onset[inc] = samples[k]
                dur[inc] = round(Int, ((samples[kn] - samples[k]) / sample_rate) * 1000)
                code[inc] = states[k]
                if haskey(d, dur[inc])
                    str = str * d[dur[inc]]
                else
                    # TODO - how to warn about this? Should ensure that ALL events have a code, or else
                    # ??? What. Should that be an error, or allowed? 
                end
                inc += 1
            end
        end
    end

    # sample onset (indices), pulse duration (in milliseconds), event code, code string
    return onset, dur, code, str
end
# ---------------------------------------------------------------------------- #
function csd_events(basedir::AbstractString)
    ifile = joinpath(basedir, "csd_events.csv")
    @assert(isfile(ifile), "Cannot locate events.csv in \"$(basedir)\"")
    ifo = CSV.parse(ifile, ',', vcat(Int, fill(Float64, 8), Int))
    good = findall(>=(0), ifo["iRew"])

    # i0 = gray->white, i1 = white->black, i2 = black->white, i3 = white->black, i4 = black->white, i5 = white->gray
    return sort!(vcat(ifo["i0"][good], ifo["i1"][good], ifo["i2"][good], ifo["i3"][good], ifo["i4"][good]))
end
# ---------------------------------------------------------------------------- #
function continuous_info(basedir::AbstractString, meta::Dict{String,Any}=metadata(basedir), nidx::Integer=typemax(Int))
    info = meta["continuous"][1]
    datadir = joinpath(basedir, "continuous", info["folder_name"])
    filepath = joinpath(datadir, "continuous.dat")
    nchan = info["num_channels"]::Int
    nsample = Int(stat(filepath).size / (nchan * sizeof(Int16)))
    fs = info["sample_rate"]::Float64

    idx = read_npy_1D(Int64, joinpath(datadir, "sample_numbers.npy"), nidx)

    return filepath, nchan, nsample, fs, idx
end
# ---------------------------------------------------------------------------- #
function read_npy_1D(::Type{T}, ifile::AbstractString, nel::Integer=typemax(Int)) where {T<:Number}
    return open(ifile, "r") do io
        magic = read(io, 6)
        @assert(String(magic) == "\x93NUMPY", "input file is not a NumPy array")

        # major = read(io, UInt8); minor = read(io, UInt8)
        seek(io, 8)

        hdr_len = read(io, UInt16)

        hdr = strip(String(read(io, hdr_len)))

        m = match(r"'shape':\s*\(\s*([^)]+)\s*\)\s*,", hdr)
        siz = parse.(Int64, split(m[1], ",", keepempty=false))

        @assert(length(siz) == 1, "NumPy array is not 1D!")

        data = Vector{T}(undef, min(siz[1], nel))

        read!(io, data)

        return data
    end
end
# ---------------------------------------------------------------------------- #
function find_channel_by_name(ifo::AbstractVector, name::AbstractString)
    for k in eachindex(ifo)
        if ifo[k]["channel_name"] == name
            return ifo[k]
        end
    end
    error("Failed to locate channel \"$(name)\"")
end
# ============================================================================ #
function find_sequence(seq::Vector{<:Integer}, target::AbstractVector{<:Integer})
    out = Vector{UnitRange{Int}}(undef, 0)
    len = length(target)
    k = 1
    while k <= (length(seq) - len) + 1
        idx = k:k+len-1
        if seq[idx] == target
            push!(out, idx)#onset[idx])
            k += len
        else
            k += 1
        end
    end
    return out
end
# ============================================================================ #
function find_sequence_re(str::String, re::Regex, len::Int)
    out=[x.offset:x.offset+len for x in eachmatch(re, str)]
    return out
end
# ============================================================================ #
end
