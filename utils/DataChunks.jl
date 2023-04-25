module DataChunks

using Mmap

export DataChunk, nchunks, ChunkDim, getchunk, swap_chunkdim
# ============================================================================ #
const TypeDict = Dict{String,DataType}("f64" => Float64, "f32" => Float32,
    "f16" => Float16, "c64" => ComplexF64, "c32" => ComplexF32,
    "c16" => ComplexF16, "i64" => Int64, "i32" => Int32, "i16" => Int16,
    "i8" => Int8, "u64" => UInt64, "u32" => UInt32, "u16" => UInt16,
    "u8" => UInt8, "b8"  => Bool)
# ============================================================================ #
@inline ChunkDim(x::Integer) = Val(x)
# ============================================================================ #
struct DataChunk{L<:AbstractArray,T,N,D} <: AbstractArray{T,N}
    x::Vector{L}
    total_size::NTuple{N,Int}
    sizes::Vector{Int}
end
# ---------------------------------------------------------------------------- #
function DataChunk(::Type{T}, total_size::NTuple{N,<:Integer}, ::Val{D}=Val(N); nchunks::Integer=Threads.nthreads()) where {T,N,D}

    x = Vector{Array{T,N}}(undef, nchunks)
    len = chunk_lengths(total_size[D], nchunks)

    for k in eachindex(x)
        siz = collect(total_size)
        siz[D] = len[k]
        x[k] = Array{T,N}(undef, siz...)
    end

    return DataChunk{Array{T,N},T,N,D}(x, total_size, cumsum(len))
end
# ---------------------------------------------------------------------------- #
function DataChunk(data::T, ::Val{D}=Val(ndims(T)); nchunks::Integer=Threads.nthreads()) where {T<:AbstractArray,D}
    siz = size(data)
    len = chunk_lengths(siz[D], nchunks)
    idx = chunk_indices(len)
    x = map(k -> copy(selectdim(data, D, k)), idx)
    L = typeof(x[1])
    return DataChunk{L, eltype(L), ndims(L), D}(x, siz, cumsum(len))
end
# ---------------------------------------------------------------------------- #
function DataChunk!(data::L, ::Val{D}=Val(ndims(L)); nchunks::Integer=Threads.nthreads()) where {L<:AbstractArray,D}
    N = ndims(L)
    siz = size(data)
    len = chunk_lengths(siz[D], nchunks)
    idx = chunk_indices(len)

    # instead of selectdim we pre-compute indices and use a view, this
    # appears to ensure type stability
    ids = map(Base.OneTo{Int}, siz)
    x = map(idx) do k
        view(data, ids[1:D-1]..., k, ids[D+1:end]...)
    end

    return DataChunk{typeof(x[1]), eltype(L), N, D}(x, siz, cumsum(len))
end
# ---------------------------------------------------------------------------- #
Base.size(d::DataChunk) = d.total_size
nchunks(d::DataChunk) = length(d.sizes)
chunks(d::DataChunk) = d.x
getchunk(d::DataChunk, k::Integer) = d.x[k]
# ---------------------------------------------------------------------------- #
function Base.getindex(d::DataChunk{L,T,N,D}, k::Vararg{Int,N}) where {L,T,N,D}
    kc = searchsortedfirst(d.sizes, k[D])
    ks = k[D] - (kc == 1 ? 0 : d.sizes[kc - 1])
    return getindex(d.x[kc], k[1:D-1]..., ks, k[D+1:end]...)
end
# ---------------------------------------------------------------------------- #
function chunk_indices(lengths::Vector{<:Integer})
    idx = Vector{UnitRange{Int}}(undef, length(lengths))
    c = 0
    for k in eachindex(lengths)
        c += lengths[k]
        n = lengths[k] - 1
        ks = (c - n)
        idx[k] = ks : ks + n
    end
    return idx
end
# ---------------------------------------------------------------------------- #
function chunk_lengths(dim_size::Integer, nchunks::Integer)
    nel = floor(Int, dim_size/nchunks)
    r = dim_size - nel*nchunks
    len = fill(nel, nchunks)
    len[1:r] .+= 1
    return len
end
# ============================================================================ #
numstr(::Type{T}) where {T<:Real} = lowercase(string(T)[1]) * string(sizeof(T)*8)
numstr(::Type{T}) where {T<:Complex} = lowercase(string(T)[1]) * string(sizeof(T)*4)
# ============================================================================ #
function format_ext(::Type{T}, ofile::AbstractString) where {T<:Number}
    m = match(r"[^\.]+((?:\.\w+)?\.\w+)", ofile)

    if m != nothing
        ofile = replace(ofile, m[1] => "." * numstr(T) * ".dat")
    else
        ofile *= "." * numstr(T) * ".dat"
    end

    return ofile
end
# ============================================================================ #
function write(ofile::AbstractString, d::DataChunk{L,T,N}) where {L,T,N}
    ofile = format_ext(T, ofile)
    open(ofile, "w") do io
        Base.write(io, UInt64(N))
        Base.write(io, map(UInt64, collect(size(d))))
        for x in chunks(d)
            Base.write(io, x)
        end
    end
    return ofile
end
# ---------------------------------------------------------------------------- #
function write(ofile::AbstractString, d::Array{T,N}) where {T,N}
    ofile = format_ext(T, ofile)
    open(ofile, "w") do io
        Base.write(io, UInt64(N))
        Base.write(io, map(UInt64, collect(size(d))))
        Base.write(io, d)
    end
    return ofile
end
# ============================================================================ #
function read_dims(ifile::AbstractString)
    return open(ifile, "r") do io
        read_dim_info(io)
    end
end
# ---------------------------------------------------------------------------- #
function read_dims(io::IOStream)
    ndim = Base.read(io, UInt64)
    siz = Vector{UInt64}(undef, ndim)
    Base.read!(io, siz)
    return siz
end
# ============================================================================ #
function read(ifile::AbstractString, ::Val{D}; nchunks::Integer=Threads.nthreads()) where D
    m = match(r"[^\.]+\.(\w+)\.dat", ifile)
    m == nothing && error("Input file extension is not valid: \"$(ifile)\"")

    T = get(TypeDict, m[1], nothing)
    T == nothing && error("Invalid file extension \"$(m[1])\"")

    return read(T, Val(D), ifile, nchunks=nchunks)
end
# ---------------------------------------------------------------------------- #
function read(::Type{T}, ::Val{D}, ifile::AbstractString; nchunks::Integer=Threads.nthreads()) where {T<:Number,D}
    !isfile(ifile) && error("Invalid input file \"$(ifile)\"")
    return open(ifile, "r") do io
        siz = read_dims(io)
        d = DataChunk(T, tuple(Int.(siz)...), Val(D), nchunks=nchunks)
        return read!(d, io, tuple(siz...))
    end
end
# ---------------------------------------------------------------------------- #
function read!(d::DataChunk, ifile::AbstractString)
    open(ifile, "r") do io
        siz = read_dims(io)
        !all(siz .== size(d)) && error("DataChunk has size $(size(d)) but file as size $(siz)")
        read!(d, io, tuple(siz...))
    end
    return d
end
# ---------------------------------------------------------------------------- #
function read!(d::DataChunk{L,T,N}, io::IOStream, siz::NTuple{N,<:Integer}) where {L,T,N}
    raw = Mmap.mmap(io, Array{T,N}, siz; grow=false, shared=false)
    k = 1
    for x in chunks(d)
        n = length(x)
        copyto!(x, 1:n, raw, k:k+n-1)
        k += n
    end
    return d
end
# ============================================================================ #
function swap_chunkdim(d::DataChunk{L,T,N,D}, ::Val{S}) where {L,T,N,D,S}
    out = DataChunk(T, size(d), Val(S), nchunks=nchunks(d))
    ax = axes(d)
    inc = 1

    for k in 1:nchunks(out)
        getchunk(out, k) .= view(d, ax[1:S-1]..., inc:out.sizes[k], ax[S+1:end]...)
        inc = out.sizes[k] + 1
    end

    return out
end
# ============================================================================ #
function pmap!(f!::Function, d::DataChunk{L,T,N,D}) where {L,T,N,D}

    @assert(nchunks(d) <= Threads.nthreads(), "Number of chunks must not be > number of threads")

    Threads.@threads for k in 1:nchunks(d)
        for slice in eachslice(getchunk(d, k), dims=D)
            f!(slice)
        end
    end

    return d
end
# ---------------------------------------------------------------------------- #
function pmap!(f!::Function, d1::DataChunk{L1,T1,N,D}, d2::DataChunk{L2,T2,N,D}) where {L1,L2,T1,T2,N,D}

    nchnk = nchunks(d1)
    @assert(nchnk <= Threads.nthreads(), "Number of chunks must not be > number of threads")
    @assert(nchnk == nchunks(d2), "Number of chunks of both data inputs must match")

    Threads.@threads for k in 1:nchnk
        for (s1,s2) in zip(eachslice(getchunk(d1,k), dims=D), eachslice(getchunk(d2,k), dims=D))
            f!(s1, s2)
        end
    end

    return d1, d2
end
# ============================================================================ #
function map!(f!::Function, d::DataChunk{L,T,N,D}) where {L,T,N,D}

    for k in 1:nchunks(d)
        for slice in eachslice(getchunk(d, k), dims=D)
            f!(slice)
        end
    end

    return d
end
# ============================================================================ #
end
