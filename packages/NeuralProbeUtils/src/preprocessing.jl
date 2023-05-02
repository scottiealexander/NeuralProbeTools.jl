# ============================================================================ #
struct Resampler <: AbstractProcessor
    rate::Rational{Int}
    fir::FIRFilter{DSP.Filters.FIRDecimator{Float64}}
end
Resampler(ratio::Rational{Int}) = Resampler(ratio, FIRFilter(DSP.resample_filter(ratio), ratio))
Base.copy(r::Resampler) = Resampler(r.rate, deepcopy(r.fir))
output_size(rr::Resampler, len::Integer, nchan::Integer) = DSP.outputlength(rr.fir, len), nchan
@inline process!(out::RealMat, inp::RealMat, r::Resampler) = resample!(out, inp, r)
# ============================================================================ #
struct Filterer <: AbstractProcessor
    coefb::Vector{Float64}
    coefa::Vector{Float64}
end
function Filterer(fs::Real, lo::Real, hi::Real)
    f = get_filter(fs, lo, hi)
    return Filterer(DSP.coefb(f), DSP.coefa(f))
end
Base.copy(f::Filterer) = Filterer(f.coefb, f.coefa)
@inline output_size(::Filterer, len::Integer, nchan::Integer) = len, nchan
@inline process!(out::RealMat, ::RealMat, f::Filterer) = filtfilt_approx!(out, f)
# ============================================================================ #
struct Interpolator <: AbstractProcessor
    bad::Vector{Int}
    n::Int
end
Base.copy(it::Interpolator) = Interpolator(it.bad, it.n)
@inline output_size(::Interpolator, len::Integer, nchan::Integer) = len, nchan
@inline process!(out::RealMat, ::RealMat, it::Interpolator) = interpolate_channels!(out, it)
# ============================================================================ #
struct DepthAverager <: AbstractProcessor
    channel_depth::Vector{Float64}
    reference_channels::Vector{Int}
end
Base.copy(d::DepthAverager) = DepthAverager(d.bad, d.n)
@inline output_size(d::DepthAverager, len::Integer, ::Integer) = len, length(unique(d.channel_positions))
@inline process!(out::RealMat, inp::RealMat, d::DepthAverager) = depth_average!(out, inp, d)
# ============================================================================ #
function get_filter(fs::Real, lo::Real, hi::Real)
    if lo > 0 && hi > 0
        f = digitalfilter(Bandpass(lo, hi, fs=fs), Butterworth(4))
    elseif lo > 0
        f = digitalfilter(Highpass(lo, fs=fs), Butterworth(4))
    elseif hi > 0
        f = digitalfilter(Lowpass(hi, fs=fs), Butterworth(4))
    else
        error("Invalid input, either <hi> or <lo> must be specified")
    end
    return f
end
# ============================================================================ #
@inline resample!(out::AbstractMatrix{Float64}, inp::AbstractMatrix, f::Resampler) = resample!(out, inp, f.fir, f.rate)

function resample!(out::AbstractMatrix{Float64}, inp::AbstractMatrix, fir::FIRFilter, rate::Rational{Int})

    @assert(size(out, 2) == size(inp, 2), "number of columns does not match")

    td = DSP.timedelay(fir)
    DSP.setphase!(fir.kernel, td)

    N = size(inp, 1)

    out_len = ceil(Int, N * rate)
    req_len = DSP.inputlength(fir, out_len)
    buf_len = DSP.outputlength(fir, req_len)

    ibuffer = zeros(Float64, req_len)
    obuffer = zeros(Float64, buf_len)

    for (k, col) in enumerate(eachcol(inp))

        # set up input buffer so it holds the data and zero-padding
        copyto!(ibuffer, 1, col, 1, N)
        ibuffer[N+1:end] .= 0.0

        if k > 1
            # reset our FIR filter state
            fill!(fir.history, 0.0)
            DSP.reset!(fir.kernel)
            DSP.setphase!(fir.kernel, td)
        end

        # apply filter and copy result (<obuffer>) to given output
        n = DSP.filt!(obuffer, fir, ibuffer)
        copyto!(view(out, :, k), 1, obuffer, 1, n)
    end

    return out
end
# ============================================================================ #
@inline filtfilt_approx!(data::AbstractMatrix{Float64}, f::Filterer) = filtfilt_approx!(data, f.coefb, f.coefa)

function filtfilt_approx!(data::AbstractMatrix{Float64}, coefb::AbstractVector{<:Real}, coefa::AbstractVector{<:Real})
    DSP.filt!(data, coefb, coefa, data)
    reverse!(data, dims=1)
    DSP.filt!(data, coefb, coefa, data)
    return reverse!(data, dims=1)
end
# ============================================================================ #
@inline interpolate_channels!(data::AbstractMatrix{<:Real}, it::Interpolator) = interpolate_channels!(data, it.bad, it.n)

function interpolate_channels!(data::AbstractMatrix{<:Real}, bad::AbstractVector{<:Integer}, n::Integer)

    nchan = size(data, 2)
    good = setdiff(1:size(data, 2), bad)

    for k = 1:length(bad)
        dist = (bad[k] .- good).^2
        ks = sortperm(dist)[1:n]

        # inverse square weighting
        tmp = (maximum(dist[ks]) + 1) .- dist[ks]
        wi = tmp ./ sum(tmp)

        data[:,bad[k]] .= (view(data, :, good[ks]) * wi)
    end

    return data
end
# ============================================================================ #
@inline depth_average!(out::AbstractMatrix{<:Real}, inp::AbstractMatrix{<:Real}, d::DepthAverager) = depth_average!(out, inp, d.channel_depth, d.reference_channels)

function depth_average!(out::AbstractMatrix{<:Real}, inp::AbstractMatrix{<:Real},
    channel_depth::AbstractVector{<:Real},
    reference_channels::AbstractVector{<:Integer}=Int[])

    ud = sort!(unique(channel_depth))

    for k in eachindex(ud)
        # find all channels with depth ud[k]
        kd = findall(isequal(ud[k]), channel_depth)

        # remove any channels w/in reference_channels
        !isempty(reference_channels) && filter!(!in(reference_channels), kd)

        # average channels w/ same depth together
        if length(kd) > 1
            out[:,k] .= mean_n(view(inp, :, kd), dim=2)
        elseif length(kd) == 1
            out[:,k] .= view(inp, :, kd[1])
        end
    end

    return out
end
# ============================================================================ #
