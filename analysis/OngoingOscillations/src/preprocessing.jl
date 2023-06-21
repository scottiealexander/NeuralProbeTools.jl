# ============================================================================ #
struct NpxProcessor <: AbstractProcessor
    resampler::Resampler
    lpfilter::Filterer
    notchfilter::Filterer
    depth_averager::DepthAverager
    buffer::Matrix{Float64}
end
# ---------------------------------------------------------------------------- #
function NpxProcessor(ratio::Rational{Int}, fs::Real, lo::Real, hi::Real, npre::Integer, npost::Integer)

    probe = Neuropixel3A()
    resampler = Resampler(ratio)
    nchan = n_channel(probe)
    len, _ = output_size(resampler, npre + npost, nchan)

    nf = iirnotch(60.0, 2.0, fs=fs)

    return NpxProcessor(resampler,
        Filterer(Int(fs * ratio), lo, hi),
        Filterer(DSP.coefb(nf), DSP.coefa(nf)),
        DepthAverager(channel_positions(probe)[:,2], reference_channels(probe)),
        Matrix{Float64}(undef, len, nchan)
    )
end
# ---------------------------------------------------------------------------- #
function Base.copy(cp::NpxProcessor)
    return NpxProcessor(copy(cp.resampler), copy(cp.filter),
        copy(cp.depth_averager), copy(cp.buffer)
    )
end
# ---------------------------------------------------------------------------- #
function output_size(p::NpxProcessor, len::Integer, nchan::Integer)
    len, nchan = output_size(p.resampler, len, nchan)
    return output_size(p.depth_averager, len, nchan)
end
# ---------------------------------------------------------------------------- #
output_type(::Type{T}, ::NpxProcessor) where T = Float64
# ---------------------------------------------------------------------------- #
function process!(out::AbstractMatrix{<:Real}, inp::AbstractMatrix{<:Real}, p::NpxProcessor)

    resample!(p.buffer, inp, p.resampler)

    median_subtract!(p.buffer, 24, p.depth_averager.reference_channels)

    filtfilt_approx!(p.buffer, p.lpfilter)
    filtfilt_approx!(p.buffer, p.notchfilter)

    depth_average!(out, p.buffer, p.depth_averager)

    return out
end
# ============================================================================ #
struct TVEProcessor <: AbstractProcessor
    pre::NpxProcessor
    freq::Vector{Float64}
    fs::Float64
    ncycle::Int
    sigma::Int
    channel_idx::Vector{Int}
    kernel::Vector{Float64}
    buffer::Matrix{Float64}
end
function TVEProcessor(pre::NpxProcessor, freq::AbstractVector{<:Real}, fs::Real, ncycle::Integer,
    sigma::Integer, channel_idx::AbstractVector{<:Integer}, window_len::Integer)

    len, nchan = output_size(p.pre, window_len, n_channel(Neuropixel3A()))

    return TVEProcessor(pre, freq, fs, ncycle, sigma, channel_idx,
        fill(1.0/window_len, window_len), Matrix{Float64}(len, nchan)
    )
end
# ---------------------------------------------------------------------------- #
Base.copy(p::TVEProcessor) = TVEProcessor(copy(p.pre), p.freq, p.fs, p.ncycle, p.sigma, p.channel_idx, p.kernel, copy(p.buffer))
output_type(::TVEProcessor) = Float64
# ---------------------------------------------------------------------------- #
function output_size(p::TVEProcessor, len::Integer, nchan::Integer)
    len, _ = output_size(p.pre, len, nchan)
    return len, 1
end
# ---------------------------------------------------------------------------- #
function process!(out::AbstractMatrix{<:Real}, inp::AbstractMatrix{<:Real}, p::TVEProcessor)
    process!(p.buffer, inp, p.pre)
    tmp = frequency_envelope!(view(p.buffer, :, p.channel_idx), p.freq, p.fs, p.ncycle. p.sigma)
    imfilter!(view(out,:,1), mean_n(tmp, dims=2), p.kernel)
    return out
end
# ============================================================================ #
