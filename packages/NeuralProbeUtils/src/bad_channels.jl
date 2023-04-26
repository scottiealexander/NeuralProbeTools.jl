# ============================================================================ #
function interpolate_bad_channels!(data::DataChunk{L,T,3,3}, bad::AbstractVector{<:Integer}, n::Integer) where {L,T}
    ri = RealInterpolator(bad, n)
    Threads.@threads for k in 1:nchunks(data)
        for slice in eachslice(getchunk(data,k), dims=3)
            interpolate_channels!(slice, ri)
        end
    end
    return data
end
# ============================================================================ #
function find_bad_channels(erp::AbstractMatrix{<:Real}, thr::Real=0.5)

    sd = apply_n(std, erp, dims=1)
    msd = maxfilter(sd, Val(5))

    bad = Int[]
    for k in eachindex(sd)
        @inbounds if ((msd[k] - sd[k]) / msd[k]) > thr
            push!(bad, k)
        end
    end
    return bad
end
# ============================================================================ #
# NOTE: **USERS MUST NOT** touch the ptr field (see cpush!() for why)
struct CircBuffer{T,N}
    x::Vector{T}
    ptr::Ref{Int}
end
CircBuffer{N}(val::T) where {N,T} = CircBuffer{T,N}(fill(val, N), Ref(1))
CircBuffer(x::Vector{T}) where {T} = CircBuffer{T,length(x)}(x, Ref(1))
# ---------------------------------------------------------------------------- #
function cpush!(c::CircBuffer{T,N}, val::T) where {T,N}
    @inbounds c.x[c.ptr[]] = val
    c.ptr[] = mod(c.ptr[], N) + 1
end
# ---------------------------------------------------------------------------- #
@generated function cmax(c::CircBuffer{T,N}) where {T,N}
    ex = Expr(:call, :max)
    for k = 1:N
        push!(ex.args, :(getindex(c.x, $k)))
    end
    return ex
end
# ============================================================================ #
function maxfilter(x::AbstractVector{T}, ::Val{N}) where {N,T<:Number}
    out = similar(x)
    @inbounds out[1] = x[1]
    buffer = CircBuffer{N}(x[1])
    for k in 2:length(out)
        cpush!(buffer, x[k])
        @inbounds out[k] = cmax(buffer)
    end
    return out
end
# ============================================================================ #
