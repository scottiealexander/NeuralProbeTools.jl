module GLX

using Mmap, NeuralProbeUtils

export spikeglx_data, event_times

# ============================================================================ #
function spikeglx_data(basedir::AbstractString, ap::Bool=false)

    apfile, lffile, metafile = all_filepaths(basedir)

    filepath = ap ? apfile : lffile
    metafile = ap ? metafile : replace(metafile, ".ap.meta"=>".lf.meta")

    meta = parse_metafile(metafile)

    if meta["typeThis"] != "imec"
        error("Invalid probe type \"$(meta["typeThis"])\"")
    end

    nchan::Int = meta["nSavedChans"]
    nsample = floor(Int, meta["fileSizeBytes"]::Float64 / (2 * nchan))
    fs::Float64 = meta["imSampRate"]

    return FlatBinaryFile{Neuropixel3A,Int16}(Neuropixel3A(), filepath, nchan, nsample, fs)
end
# ============================================================================ #
struct GLXData
    apfile::String
    lffile::String
    ap_meta::Dict{String, Any}
    lf_meta::Dict{String, Any}
end
# ============================================================================ #
function GLXData(ifile::AbstractString)
    apfile, lffile, metafile = all_filepaths(ifile)

    ap_meta = parse_metafile(metafile)

    if ap_meta["typeThis"] != "imec"
        error("Invalid probe type \"$(ap_meta["typeThis"])\"")
    end

    lf_meta = parse_metafile(replace(metafile, ".ap.meta"=>".lf.meta"))

    return GLXData(apfile, lffile, ap_meta, lf_meta)
end
# ============================================================================ #

ap_channels(g::GLXData) = Int(g.ap_meta["nSavedChans"])::Int - 1
lf_channels(g::GLXData) = Int(g.lf_meta["nSavedChans"])::Int - 1
ap_samples(g::GLXData) = floor(Int, g.ap_meta["fileSizeBytes"]::Float64 / (2 * (ap_channels(g) + 1)))
lf_samples(g::GLXData) = floor(Int, g.lf_meta["fileSizeBytes"]::Float64 / (2 * (lf_channels(g) + 1)))

ap_sample_rate(g::GLXData) = g.ap_meta["imSampRate"]::Float64
lf_sample_rate(g::GLXData) = g.lf_meta["imSampRate"]::Float64

# same for ap and lf
volt_scale(g::GLXData) = g.ap_meta["imAiRangeMax"]::Float64 / 512

function ap_memmap(g::GLXData)
    return open(g.apfile, "r") do io
        Mmap.mmap(io, Matrix{Int16}, (ap_channels(g) + 1, ap_samples(g)); grow=false)
    end
end
function lf_memmap(g::GLXData)
    return open(g.lffile, "r") do io
        Mmap.mmap(io, Matrix{Int16}, (lf_channels(g) + 1, lf_samples(g)); grow=false)
    end
end
# ============================================================================ #
function original_channels(g::GLXData)

    items = vcat(
        split(g.ap_meta["snsSaveChanSubset"], ','),
        split(g.lf_meta["snsSaveChanSubset"], ',')
    )

    d = Vector{Int}()
    for item in items
        tmp = tryparse(Int, item)
        if tmp != nothing
            push!(d, tmp)
        else
            m = match(r"(\d+)\:(\d+)", item)
            if m != nothing
                try
                    rn = tryparse(Int, m[1]):tryparse(Int, m[2])
                    append!(d, rn)
                catch err
                    @warn(err)
                end
            else
                error("Invalid array entry: \"$(item)\"")
            end
        end
    end

    return d .+ 1
end
# ============================================================================ #
"ap, lf, sy = channel_counts(g)"
function channel_counts(g::GLXData)
    # for ap channel
    m = tryparse.(Int, split(g.ap_meta["snsApLfSy"], ','))

    # add counts for lf channel
    m .+= tryparse.(Int, split(g.lf_meta["snsApLfSy"], ','))

    #      ap    lf    sy
    return m[1], m[2], m[3]
end
# ============================================================================ #
"ap_gain, lf_gain = channel_gain(::Type{T<:Real}, g)"
function channel_gain(::Type{T}, g::GLXData) where T<:Real
    done = false
    k = 1
    gains = Vector{Vector{T}}()
    while !done
        # imroTbl is the same regardless of ap or lf file
        m = match(r"\(([^\)]+)\)", g.ap_meta["imroTbl"], k)
        if m == nothing
            done = true
        else
            # do nothing with header information? not sure what this is for...
            if k != 1
                push!(gains, tryparse.(T, split(m[1], [' ', ','])))
            end
            k = m.offsets[1] + length(m[1]) + 1
        end
    end

    ap_gain = zeros(T, length(gains))
    lf_gain = zeros(T, length(gains))
    for k in eachindex(gains)
        ap_gain[k] = gains[k][4]
        lf_gain[k] = gains[k][5]
    end

    return ap_gain, lf_gain
end
# ============================================================================ #
function shank_map(g::GLXData)
    k = 1
    col = Vector{Int}()
    row = Vector{Int}()
    done = false
    while !done
        # snsShankMap appears to only exist in ap file metadata
        m = match(r"\(([^\)]+)\)", g.ap_meta["snsShankMap"], k)
        if m == nothing
            done = true
        else
            # do nothing with header information? not sure what this is for...
            if k != 1
                tmp = tryparse.(Float64, split(m[1], [' ', ',', ':']))

                # col / row indicators are 0 based in glx files
                push!(col, tmp[2] + 1)
                push!(row, tmp[3] + 1)
            end
            k = m.offsets[1] + length(m[1]) + 1
        end
    end
    return row, col
end
# ============================================================================ #
"""
ap_gain_correct!(data::Matrix{<:AbstractFloat}, g::GLXData)
"""
function ap_gain_correct!(data::Matrix{T}, g::GLXData) where T<:AbstractFloat
    ap_gain, _ = channel_gain(T, g)
    (size(data, 1) - 1) != length(ap_gain) && error("Data matrix has incorrect number of channels")
    return gain_correct!(data, g, ap_gain)
end
"""
lf_gain_correct!(data::Matrix{<:AbstractFloat}, g::GLXData)
"""
function lf_gain_correct!(data::Matrix{T}, g::GLXData) where T<:AbstractFloat
    _, lf_gain = channel_gain(T, g)
    (size(data, 1) - 1) != length(lf_gain) && error("Data matrix has incorrect number of channels")
    return gain_correct!(data, g, lf_gain)
end
function gain_correct!(data::Matrix{T}, g::GLXData, gains::Vector{T}) where T<:AbstractFloat
    i2v = volt_scale(g)

    # do not correct the event channel (last row)
    for k in 1:(size(data, 1) - 1)
        data[k,:] .*= (i2v / gains[k])
    end
    return data
end
# ============================================================================ #
"""
gain_correct!(data::Vector{Float64}, g::GLXData, idx::Integer)
"""
function gain_correct!(data::Vector{T}, g::GLXData, idx::Integer) where T<:AbstractFloat
    chans = original_channels(g)
    ap_gain, lf_gain = channel_gain(T, g)
    nap = length(ap_gain)
    nu = nap * 2

    i2v = volt_scale(g)

    j = Int(chans[idx])
    if j <= nap
        conv = i2v / ap_gain[j]
    elseif j <= nu
        conv = i2v / lf_gain[j - nap]
    else
        # probably the Sy channel
        conv = 1.0
    end

    data .*= conv

    return data
end
# ============================================================================ #
function all_filepaths(ifile::AbstractString)
    if isdir(ifile)
        dir, name = splitdir(ifile)
        ifile = joinpath(dir, name, name * "_t0.imec.ap.bin")
    end
    m = match(r".*\.imec.(\w+)\.(\w+)$", ifile)
    m == nothing && error("Failed to parse file name!")
    if m[1] == "ap"
        if m[2] == "bin"
            apfile = ifile
            metafile = replace(ifile, ".bin"=>".meta")
            lffile = replace(ifile, ".ap."=>".lf.")
        elseif m[2] == "meta"
            apfile = replace(ifile, ".meta"=>".bin")
            lffile = replace(ifile, ".ap."=>".lf.")
            metafile = ifile
        end

    elseif m[1] == "lf"
        if m[2] == "bin"
            lffile = ifile
            metafile = replace(ifile, ".lf.bin"=>".ap.meta")
            apfile = replace(ifile, ".lf."=>".ap.")
        elseif m[2] == "meta"
            lffile = replace(ifile, ".meta"=>".bin")
            apfile = replace(ifile, ".lf."=>".ap.")
            metafile = ifile
        end
    end
    return apfile, lffile, metafile
end
# ============================================================================ #
function parse_metafile(ifile::AbstractString)
    d = Dict{String, Any}()
    for line in eachline(ifile)
        m = match(r"([^\=]+)\s*\=\s*([^\n]*)", line)
        if m != nothing
            k = strip(m[1], ['\t','\n','\v','\f','\r',' ','~'])
            tmp = strip(m[2])
            v = tryparse(Float64, tmp)
            if v != nothing
                d[k] = v
            else
                d[k] = tmp
            end
        end
    end

    return d
end
# ============================================================================ #
"""
data = read_ap_data(::Type{<:AbstractFloat}, g::GLXData, start::Integer=0, length::Integer=-1)
"""
function read_ap_data(::Type{T}, g::GLXData, start::Integer=0, length::Integer=-1) where T<:AbstractFloat
    return ap_gain_correct!(
        _read_data(T, g.apfile, ap_channels(g) + 1, ap_samples(g), start, length),
        g
    )
end
"""
data = read_ap_data(g::GLXData, start::Integer=0, length::Integer=-1)
"""
function read_ap_data(g::GLXData, start::Integer=0, length::Integer=-1)
    return _read_data(Int16, g.apfile, ap_channels(g) + 1, ap_samples(g),
        start, length)
end
# ============================================================================ #
"""
data = read_lf_data(::Type{<:AbstractFloat}, g::GLXData, start::Integer=0, length::Integer=-1)
"""
function read_lf_data(::Type{T}, g::GLXData, start::Integer=0, length::Integer=-1) where T<:AbstractFloat
    return lf_gain_correct!(
        _read_data(T, g.lffile, lf_channels(g) + 1, lf_samples(g), start, length),
        g
    )
end
"""
data = read_lf_data(g::GLXData, start::Integer=0, length::Integer=-1)
"""
function read_lf_data(g::GLXData, start::Integer=0, length::Integer=-1)
    return _read_data(Int16, g.lffile, lf_channels(g) + 1, lf_samples(g),
        start, length)
end
# ============================================================================ #
function _read_data(::Type{T}, ifile::AbstractString, nchan::Integer,
    nsamp::Integer, start::Integer, length::Integer) where T<:Real

    length = length < 0 ? nsamp : length
    nsamp = min(length, nsamp - start)

    offset = Int(max(0, start) * sizeof(Int16) * nchan)

    data = Matrix{Int16}(undef, Int(nchan), Int(nsamp))

    open(ifile, "r") do io
        seek(io, offset)
        read!(io, data)
    end

    return Matrix{T}(data)
end
# ============================================================================ #
"""
data = read_channel(g::GLXData, idx::Integer)
"""
function read_channel(g::GLXData, idx::Integer)::Vector{Int16}
    event = false
    chans = original_channels(g)

    kchan = chans[idx]

    # chans[385] and chans[770] both == 769 (event channel)
    if kchan > (ap_channels(g) + lf_channels(g))
        # event channel
        kchan = lf_channels(g) + 1
        return lf_memmap(g)[kchan,:]
    elseif kchan > ap_channels(g)
        # lf channel
        kchan -= ap_channels(g)
        return lf_memmap(g)[kchan,:]
    else
        # ap channel
        return ap_memmap(g)[kchan,:]
    end
end
# ============================================================================ #
"""
***NOTE***: if a sy.bin file already exists this will do nothing and return
an empty Int16 vector
"""
function write_sync_file(g::GLXData)
    ofile = replace(g.apfile, "imec.ap.bin" => "imec.sy.bin")

    # Windows and Mmap do not play nicely together, so if a .sy.bin file already
    # exists just given a warning and punt
    if !isfile(ofile)
        raw = ap_memmap(g)
        d = raw[end,:] .- raw[end,1]
        out = Mmap.mmap(ofile, Vector{Int16}, (ap_samples(g),); grow=false, shared=false)
        out .= d
        Mmap.sync!(out)
        return d
    else
        println("[WARNING]: File \"", ofile, "\" already exists, remove manually to re-create it")
        return Int16[]
    end
end
# ============================================================================ #
function event_channel(g::GLXData)
    ifile = replace(g.apfile, "imec.ap.bin" => "imec.sy.bin")
    if isfile(ifile)
        return open(ifile, "r") do io
            Mmap.mmap(io, Vector{Int16}, (ap_samples(g),); grow=false, shared=false)
        end
    else
        return write_sync_file(g)
    end
end
# ============================================================================ #
function event_indices(evt::Vector{Int16}, bit::Integer)
    !(1 <= bit <= 16) && error("Invalid bit: must be between 1 and 16")
    onset = Vector{Int}()
    offset = Vector{Int}()
    last_val = bitget(evt[1], bit)
    @inbounds for k in 2:length(evt)
        cur_val = bitget(evt[k], bit)
        df = cur_val - last_val
        # if  cur_val - last_val > 0
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
# ============================================================================ #
function event_times(g::GLXData, bit::Integer)
    sps = 1.0 / ap_sample_rate(g)
    onset, _ = event_indices(event_channel(g), bit)
    return (onset .- 1) .* sps
end
# ============================================================================ #
@inline bitget(x::Integer, k::Integer) = (x >> (k-1)) & 0x01
# ============================================================================ #
function write_data(g::GLXData, data::AbstractArray{<:Real,3},
    onset::AbstractVector{<:Integer}; ofile::AbstractString="",
    suffix::AbstractString="", ap::Bool=true, force::Bool=false)

    glxfile = ap ? g.apfile : g.lffile
    type_str = ap ? "ap" : "lf"
    if isempty(ofile)
        odir = joinpath(dirname(glxfile), "artefact_correct")
        !isdir(odir) && mkdir(odir)
        ofile = joinpath(
            odir,
            replace(basename(glxfile), ".imec."*type_str*".bin" => suffix * ".imec."*type_str*".bin")
        )
    elseif (ofile == g.lffile || ofile == g.apfile) && !force
        error("Output and input files are the same, use option \"force=true\" to force output")
    end

    #@info("ofile = \"$(ofile)\"")

    # delete output file if it already exists to avoid clashes
    if isfile(ofile)
        #@info("Removing existing output file")
        # sometimes Mmap files don't get closed (causing an error in rm() /
        # unlink()), so if that happens try opening then closing the file
        try
            rm(ofile, force=true)
        catch err
            io = open(ofile, "r")
            close(io)
            rm(ofile, force=true)
        end
    end

    # first copy the file from which our tensor came (`data`) and then...
    #@info("Copying " * type_str * "file...")
    if !copy_file(glxfile, ofile)
        error("Failed to copy file!")
    end

    nchan = ap ? ap_channels(g) : lf_channels(g)
    nsamp = ap ? ap_samples(g) : lf_samples(g)

    # overwrite the artefact contaminated segments "in place"
    open(ofile, "r+") do io
        mm = Mmap.mmap(io, Matrix{Int16}, (nchan + 1, nsamp); grow=false, shared=false)

        len = size(data, 1) - 1
        kchan = 1:size(data, 2)

        for k in eachindex(onset)
            idx = onset[k]:onset[k] + len
            mm[kchan, idx] .= round.(Int16, data[:,:,k]')
        end

        Mmap.sync!(mm)
    end

    # file extension regex pattern that covers both ap and lf cases
    ext_ptrn = r".imec.(?:lf|ap).bin"

    # NOTE: we no longer copy the lf data when the ap data is corrected, as that
    # needs to be corrected seperatly and having an "_ac" version will cause
    # problems

    # copy sync file, note that we use copy_file() instead of safe_cp() as
    # cp() fails with large (or binary?) files on windows
    syncfile = replace(ofile, ext_ptrn => ".imec.sy.bin")
    !isfile(syncfile) && copy_file(replace(g.apfile, ".ap." => ".sy."), syncfile)

    # copy the .meta file for the ap data
    safe_cp(replace(g.apfile, ".bin" => ".meta"), replace(ofile, ext_ptrn => ".imec.ap.meta"), force=false)

    # copy the .meta file for the lf data
    safe_cp(replace(g.lffile, ".bin" => ".meta"), replace(ofile, ext_ptrn => ".imec.lf.meta"), force=false)

    return ofile
end
# ============================================================================ #
function copy_file(src::AbstractString, dst::AbstractString)
    @static if Sys.iswindows()
        run(`cmd /c copy /B $(src) /Y $(dst) ">" nul`)
    else
        cp(src, dst, force=true)
    end
    return isfile(dst)
end
# ============================================================================ #
function safe_cp(src::AbstractString, dst::AbstractString; force::Bool=false)
    if !isfile(dst)
        cp(src, dst, force=false)
    elseif force
        cp(src, dst, force=true)
    end
    return nothing
end
# ============================================================================ #
end
