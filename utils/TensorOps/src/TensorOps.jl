module TensorOps

using LinearAlgebra, Statistics, ImageFiltering

include("./csd.jl")

export mean_3, std_3, mean_n, apply_n, rm_baseline!, median_subtract!, smooth, smooth!
# ============================================================================ #
@inline apply_n(f::Function, x::AbstractArray; dims::Integer) = dropdims(f(x, dims=dims), dims=dims)
mean_n(x::AbstractArray; dims::Integer) = apply_n(mean, x, dims=dims)

mean_3(x::AbstractArray{T,3}) where {T} = apply_n(mean, x, dims=3)
std_3(x::AbstractArray{T,3}) where {T} = apply_n(std, x, dims=3)
# ============================================================================ #
function demean!(x::AbstractMatrix{<:Real})
    for col in eachcol(x)
        col .-= mean(col)
    end
    return x
end
# ============================================================================ #
function rm_baseline!(x::AbstractMatrix{<:Real}, nbl::Integer)
    for col in eachcol(x)
        col .-= mean(view(col, 1:nbl))
    end
    return x
end

function rm_baseline!(x::AbstractArray{<:Real, 3}, nbl::Integer)
    for slice in eachslice(x, dims=3)
        rm_baseline!(slice, nbl)
    end
    return x
end
# ============================================================================ #
function median_subtract!(x::AbstractMatrix{<:Real})
    x .-= median(x, dims=2)
    return x
end
function median_subtract!(x::AbstractArray{<:Real,3})
    for slice in eachslice(x, dims=3)
        median_subtract!(slice)
    end
    return x
end
# for Neuropixel phase3A, step = 24, omit = REFERENCE_CHANNELS
function median_subtract!(x::AbstractMatrix{<:Real}, step::Integer, omit::AbstractVector{<:Integer}=Int[])
    for k in 1:step
        kchan = setdiff(k:step:size(x, 2), omit)
        md = median(view(x, :, kchan), dims=2)
        x[:,kchan] .-= md
    end
    return x
end
function median_subtract!(x::AbstractArray{<:Real,3}, step::Integer, omit::AbstractVector{<:Integer}=Int[])
    for slice in eachslice(x, dims=3)
        median_subtract!(slice, step, omit)
    end
    return x
end
# ============================================================================ #
const KT = KernelFactors.ReshapedOneD
smooth(x::AbstractMatrix, kernels::Tuple{KT,KT}) = imfilter(x, kernels)
smooth!(x::AbstractMatrix, kernels::Tuple{KT,KT}) = imfilter!(x, x, kernels)

smooth(x::AbstractMatrix, sd::Tuple{Int,Int}) = smooth(x, KernelFactors.IIRGaussian(sd))
smooth(x::AbstractMatrix, sd::Vector{<:Integer}) = smooth(x, KernelFactors.IIRGaussian((sd[1], sd[2])))
smooth!(x::AbstractMatrix, sd::Tuple{Int,Int}) = smooth!(x, KernelFactors.IIRGaussian(sd))
smooth!(x::AbstractMatrix, sd::Vector{<:Integer}) = smooth!(x, KernelFactors.IIRGaussian((sd[1], sd[2])))

function smooth!(x::AbstractArray{<:Real,3}, sd::Tuple{Int,Int})
    kernels = KernelFactors.IIRGaussian(sd)
    for slice in eachsilce(x, dims=3)
        smooth!(slice, kernels)
    end
    return x
end
# ============================================================================ #
end
