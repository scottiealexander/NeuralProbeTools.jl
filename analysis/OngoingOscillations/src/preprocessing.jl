# ============================================================================ #
struct NpxProcessor <: AbstractProcessor
    resampler::Resampler
    filter::Filterer
    depth_averager::DepthAverager
    buffer::Matrix{Float64}
end
# ============================================================================ #
function NpxProcessor(ratio::Rational{Int}, fs::Real, lo::Real, hi::Real, npre::Integer, npost::Integer)

    probe = Neuropixel3A()
    resampler = Resampler(ratio)
    nchan = n_channel(probe)
    len, _ = output_size(resampler, npre + npost, nchan)

    return NpxProcessor(
        resampler,
        Filterer(Int(fs * ratio), lo, hi),
        DepthAverager(channel_positions(probe)[:,2], reference_channels(probe)),
        Matrix{Float64}(undef, len, nchan)
    )
end
# ============================================================================ #
function Base.copy(cp::NpxProcessor)
    return NpxProcessor(copy(cp.resampler), copy(cp.filter),
        copy(cp.depth_averager), copy(cp.buffer)
    )
end
# ============================================================================ #
function output_size(p::NpxProcessor, len::Integer, nchan::Integer)
    len, nchan = output_size(p.resampler, len, nchan)
    return output_size(p.depth_averager, len, nchan)
end
# ============================================================================ #
output_type(::Type{T}, ::NpxProcessor) = Float64
# ============================================================================ #
function process!(out::AbstractMatrix{<:Real}, inp::AbstractMatrix{<:Real}, p::NpxProcessor)

    resample!(p.buffer, inp, p.resampler)

    median_subtract!(p.buffer, 24, p.depth_averager.reference_channels)

    filtfilt_approx!(p.buffer, p.filterer)

    depth_average!(out, p.buffer, p.depth_averager)

    return out
end
# ============================================================================ #
