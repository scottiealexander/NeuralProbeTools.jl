module Intan

using StringEncodings, Mmap

using NeuralProbeUtils

# export SingleFile, OneFilePerSignalType, OneFilePerChannel
export preprocessor, load_and_process, channel_order

const TimeType = Int32
const RawType = Int16
const IOType = UInt16

# vertical contact spacing is 25 um
const CONTACT_SPACING = 0.000025

include("./events.jl")

# ============================================================================ #
function intan_data(filepath::AbstractString, p::NeuralProbeUtils.AbstractProbe=DBCDeepArray())
    nc, ns = get_datasize(filepath)
    fs = read_sample_rate(filepath)
    return FlatBinaryFile{typeof(p),Int16}(p, filepath, nc, ns, fs)
end
# ============================================================================ #
function read_channel_positions(filepath::AbstractString)

    if isdir(filepath)
        ifile = joinpath(filepath, "channel_positions.f64.dat")
    elseif !endswith(filepath, "channel_positions.f64.dat")
        ifile = joinpath(dirname(filepath), "channel_positions.f64.dat")
    else
        ifile = filepath
    end

    @assert(isfile(ifile), "failed to locate channel_positions file!")

    return open(ifile, "r") do io
        ndim = read(io, UInt64)
        siz = Vector{UInt64}(undef, ndim)
        read!(io, siz)
        data = Array{Float64,Int(ndim)}(undef, Int.(siz)...)
        read!(io, data)
    end
end

# de-interleave channel ordering
function get_channel_order(idir::AbstractString)
    pos = read_channel_positions(idir)
    return sortperm(pos[:,2])
end

reorder!(data::AbstractVector{<:Real}, order::AbstractVector{<:Integer}) = copy!(data, data[order])
reorder!(data::AbstractMatrix{<:Real}, order::AbstractVector{<:Integer}) = copy!(data, data[:,order])

function reorder!(data::Array{T,3}, order::AbstractVector{<:Integer}) where T<:Real
    tmp = zeros(T, size(data, 1), size(data, 2))
    for slice in eachslice(data, dims=3)
        copy!(tmp, view(slice, :, order))
        copy!(slice, tmp)
    end
    return data
end

function read_time_file(filepath::AbstractString, fs::Real=-1)

    fs = read_sample_rate(filepath, fs)

    nsample = Int(stat(filepath).size / sizeof(TimeType))
    data = Vector{TimeType}(undef, nsample)
    read!(filepath, data)
    return data ./ fs
end

function read_data_file(filepath::AbstractString, nchannel::Integer=-1, nsample::Integer=-1)
    if nchannel < 1
        hdr = read_rhd_header(get_rhd_file(filepath), OneFilePerSignalType)
        nchannel = hdr["signal_group"]["n_amp_channel"]
    end

    if nsample < 1
        total_samples = Int(stat(filepath).size / sizeof(RawType))
        nsample = Int(total_samples / nchannel)
    end

    data = Matrix{RawType}(undef, nchannel, nsample)

    read!(filepath, data)

    return data
end

function get_datasize(filepath::AbstractString)
    hdr = read_rhd_header(get_rhd_file(filepath), OneFilePerSignalType)
    nchannel = Int(hdr["signal_group"]["n_amp_channel"])

    total_samples = Int(stat(filepath).size / sizeof(RawType))
    nsample = Int(total_samples / nchannel)

    return nchannel, nsample
end

function read_digital_file(filepath::AbstractString)
    nsample = Int(stat(filepath).size / sizeof(IOType))
    data = Vector{IOType}(undef, nsample)
    read!(filepath, data)

    return data
end

function event_indices(evt::Vector{IOType}, bit::Integer)
    !(1 <= bit <= 16) && error("Invalid bit: must be between 1 and 16")
    onset = Vector{Int}()
    offset = Vector{Int}()
    last_val = Int(bitget(evt[1], bit))
    @inbounds for k in 2:length(evt)
        cur_val = Int(bitget(evt[k], bit))
        df = cur_val - last_val
        if df > 0
            # rising edge
            push!(onset, k)
        elseif df < 0
            # falling edge
            push!(offset, k)
        end
        last_val = cur_val
    end
    return onset, offset
end

function event_times(evt::Vector{IOType}, bit::Integer, fs::Real)
    sps = 1.0 / fs
    onset, _ = event_indices(evt, bit)
    return (onset .- 1) .* sps
end

function event_times(evt::Vector{IOType}, t::Vector{TimeType}, bit::Integer)
    onset, _ = event_indices(evt, bit)
    return t[onset]
end

@inline bitget(x::Integer, k::Integer) = (x >> (k-1)) & 0x01

@inline get_rhd_file(filepath::AbstractString) = joinpath(dirname(filepath), "info.rhd")

function read_sample_rate(filepath::AbstractString, fs::Real=0)
    if fs <= 0
        rhd_file = get_rhd_file(filepath)
        if !isfile(rhd_file)
            @warn("Failed to locate \"info.rhd\" and sample rate was not given")
            fs = 1.0f0
        else
            fs = open(rhd_file, "r") do io
                read_sample_rate(io)
            end
        end
    else
        fs = Float32(fs)
    end
    return fs
end

function read_sample_rate(io::IOStream)
    pos = position(io)
    seek(io, 0)
    assert_intan(io)
    seek(io, sizeof(UInt32) + 2 * sizeof(Int16))
    fs = read(io, Float32)
    seek(io, pos)
    return fs
end

@inline function assert_intan(io::IOStream)
    # make sure it's really and intan rhd header
    magic_number = read(io, UInt32)
    @assert(magic_number == 0xC6912702, "file stream appears not to contain an Intan header")
    return nothing
end

@enum DataFormat SingleFile OneFilePerSignalType OneFilePerChannel
function read_rhd_header(filepath::AbstractString, mode::DataFormat=OneFilePerSignalType)
    return open(filepath, "r") do io

        hdr = Dict{String,Any}()

        assert_intan(io)

        hdr["version"] = VersionNumber(read(io, Int16), read(io, Int16))

        hdr["sample_rate"] = read(io, Float32)
        hdr["offset_filter_enabled"] = read(io, Int16)

        hdr["offset_filter_cutoff"] = read(io, Float32)
        hdr["hp_cutoff"], hdr["lp_cutoff"] = read(io, Float32), read(io, Float32)

        hdr["req_offset_cutoff"], hdr["req_hp_cutoff"], hdr["req_lp_cutoff"] = read(io, Float32), read(io, Float32), read(io, Float32)

        hdr["notch_filter"] = read(io, Int16)

        hdr["req_imped"], hdr["act_imped"] = read(io, Float32), read(io, Float32)

        hdr["notes"] = join([read_qstring(io), read_qstring(io), read_qstring(io)], '\n')

        hdr["n_temp_sensor"], hdr["board_mode"], hdr["reference_channel"] = Int16(0), Int16(0), ""

        if hdr["version"] >= v"1.1"
            hdr["n_temp_sensor"] = read(io, Int16)
        end

        if hdr["version"] >= v"1.3"
            hdr["board_mode"] = read(io, Int16)
        end

        if hdr["version"] >= v"2.0"
            hdr["reference_channel"] = read_qstring(io)
        end

        n_signal_grp = read(io, Int16)

        if mode == SingleFile

            hdr["signal_groups"] = Vector{Dict{String,Any}}(undef, n_signal_grp)

            for k = 1:n_signal_grp
                hdr["signal_groups"][k] = read_signal_group(io, true)
            end

        elseif mode == OneFilePerSignalType
            hdr["signal_group"] = read_signal_group(io, false)

        else
            error("OneFilePerChannel mode is not yet implemented...")
        end

        return hdr
    end
end

function read_signal_group(io::IOStream, read_list::Bool)
    grp = Dict{String,Any}()
    grp["name"] = read_qstring(io)
    grp["prefix"] = read_qstring(io)
    grp["enabled"] = read(io, Int16)
    grp["n_channel"] = read(io, Int16)
    grp["n_amp_channel"] = read(io, Int16)
    if read_list && grp["enabled"] > 0 && grp["n_channel"] > 0
        grp["channels"] = read_channel_list(io, grp["n_channel"])
    end
    return grp
end

function read_channel_list(io::IOStream, nchannel::Integer)
    chans = Vector{Dict{String,Any}}(undef, nchannel)
    for k in 1:nchannel
        chans[k] = Dict{String, Any}()
        chans[k]["name"] = read_qstring(io)
        chans[k]["custom_name"] = read_qstring(io)

        chans[k]["order"] = read(io, Int16)
        chans[k]["custom_order"] = read(io, Int16)
        chans[k]["signal_type"] = read(io, Int16)
        chans[k]["enabled"] = read(io, Int16)
        chans[k]["chip_channel"] = read(io, Int16)
        chans[k]["board_stream"] = read(io, Int16)
        chans[k]["trigger_mode"] = read(io, Int16)
        chans[k]["threshold"] = read(io, Int16)
        chans[k]["polarity"] = read(io, Int16)
        chans[k]["impedance"] = read(io, Float32)
        chans[k]["impedance_phase"] = read(io, Float32)
    end

    return chans
end

# read a Qt formated "q-string"
function read_qstring(io::IOStream)
    nbytes = read(io, UInt32)
    if nbytes == 0xFFFFFFFF
        return ""
    else
        return decode(read(io, nbytes), "UTF-16")
    end
end

end # module
