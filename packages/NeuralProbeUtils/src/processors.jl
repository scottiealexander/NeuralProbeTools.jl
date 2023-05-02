# ============================================================================ #
abstract type AbstractProcessor end
Base.copy(::AbstractProcessor) = error("Not implemented")
output_size(::AbstractProcessor, ::Integer, ::Integer) = error("Not implemented")
output_type(::Type{T}, ::AbstractProcessor) where T = T
process!(::RealMat, ::RealMat, ::AbstractProcessor) = error("Not implemented")
# ============================================================================ #
struct NullProcessor <: AbstractProcessor end
Base.copy(n::NullProcessor) = n
output_size(::NullProcessor, len::Integer, nchan::Integer) = len, nchan
@inline process!(out::RealMat, inp::RealMat, ::NullProcessor) = out
# ============================================================================ #
struct Preprocessor{N,T} <: AbstractProcessor
    procs::T
end
Preprocessor(procs::Vararg{<:AbstractProcessor}) = Preprocessor{length(procs),typeof(procs)}(procs)
function output_size(g::Preprocessor{N}, len::Integer, nchan::Integer) where N
    for k = 1:N
        len, nchan = output_size(g.procs[k], len, nchan)
    end
    return len, nchan
end
output_type(::Type{T}, g::Preprocessor) where T = Float64
@generated function Base.copy(g::Preprocessor{N}) where N
    ex = Expr(:call, Preprocessor)
    for k in 1:N
        push!(ex.args, :(copy(getindex(g.procs, $k))))
    end

    return ex
end
@generated function process!(out::RealMat, inp::RealMat, g::Preprocessor{N}) where N
    ex = Expr(:block)
    for k in 1:N
        push!(ex.args, :(process!(out, inp, getindex(g.procs, $k))))
    end

    return ex
end
# ============================================================================ #
