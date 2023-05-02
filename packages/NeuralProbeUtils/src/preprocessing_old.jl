# ============================================================================ #
abstract type AbstractResampler end

struct NullResampler <: AbstractResampler end

struct RealResampler <: AbstractResampler
    rate::Rational{Int}
    fir::FIRFilter{DSP.Filters.FIRDecimator{Float64}}
end

RealResampler(ratio::Rational{Int}) = RealResampler(ratio, FIRFilter(DSP.resample_filter(ratio), ratio))
output_length(::NullResampler, inp_len::Integer) = inp_len
output_length(rr::RealResampler, inp_len::Integer) = DSP.outputlength(rr.fir, inp_len)
# ============================================================================ #
abstract type AbstractFilter end

struct NullFilter <: AbstractFilter end

struct RealFilter <: AbstractFilter
    coefb::Vector{Float64}
    coefa::Vector{Float64}
end

function RealFilter(fs::Real, lo::Real, hi::Real)
    f = get_filter(fs, lo, hi)
    return RealFilter(DSP.coefb(f), DSP.coefa(f))
end
# ============================================================================ #
abstract type AbstractInterpolator end

struct NullInterpolator <: AbstractInterpolator end

struct RealInterpolator <: AbstractInterpolator
    bad::Vector{Int}
    n::Int
end
# ============================================================================ #
struct Preprocessor{R<:AbstractResampler,F<:AbstractFilter,I<:AbstractInterpolator}
    resampler::R
    filter::F
    interpolator::I
end

function Preprocessor(;
    resampler::AbstractResampler=NullResampler(),
    filter::AbstractFilter=NullFilter(),
    interpolator::AbstractInterpolator=NullInterpolator())
    return Preprocessor(resampler, filter, interpolator)
end
output_length(p::Preprocessor, inp_len::Integer) = output_length(p.resampler, inp_len)
output_type(::Type{T}, p::Preprocessor{NullResampler,NullFilter,NullInterpolator}) where {T} = T
output_type(::Type{T}, p::Preprocessor) where {T} = Float64
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
# no preprocessing
preprocessor() = Preprocessor()

# filter only
preprocessor(fs::Real, lo::Real, hi::Real) = Preprocessor(filter=RealFilter(fs, lo, hi))

# resample only
preprocessor(ratio::Rational{Int}) = Preprocessor(resampler=RealResampler(ratio))

# interpolate only
preprocessor(bad::Vector{Int}, n::Integer) = Preprocessor(interpolator=RealInterpolator(bad, n))

# resample and filter
preprocessor(ratio::Rational{Int}, fs::Real, lo::Real, hi::Real) = Preprocessor(resampler=RealResampler(ratio), filter=RealFilter(fs, lo, hi))

# resample and interpolate
preprocessor(ratio::Rational{Int}, bad::Vector{Int}, n::Integer) = Preprocessor(resampler=RealResampler(ratio), interpolator=RealInterpolator(bad, n))

# filter and interpolate
preprocessor(fs::Real, lo::Real, hi::Real, bad::Vector{Int}, n::Integer) =
    Preprocessor(filter=RealFilter(fs, lo, hi), interpolator=RealInterpolator(bad, n))

# all three
function preprocessor(ratio::Rational{Int}, fs::Real, lo::Real, hi::Real, bad::Vector{Int}, n::Integer)
    return Preprocessor(
        resampler=RealResampler(ratio),
        filter=RealFilter(fs, lo, hi),
        interpolator=RealInterpolator(bad, n)
    )
end
# ============================================================================ #
function preprocess!(out::AbstractMatrix{Float64}, inp::AbstractMatrix, proc::Preprocessor)
    resample!(out, inp, proc.resampler)
    filtfilt_approx!(out, proc.filter)
    interpolate_channels!(out, proc.interpolator)
    return out
end
# ============================================================================ #
@inline resample!(out::AbstractMatrix{Float64}, inp::AbstractMatrix, ::NullResampler) = copy!(out, inp)

resample!(out::AbstractMatrix{Float64}, inp::AbstractMatrix, f::RealResampler) = resample!(out, inp, f.fir, f.rate)

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
@inline filtfilt_approx!(data::AbstractMatrix{Float64}, ::NullFilter) = data

filtfilt_approx!(data::AbstractMatrix{Float64}, f::RealFilter) = filtfilt_approx!(data, f.coefb, f.coefa)

function filtfilt_approx!(data::AbstractMatrix{Float64}, coefb::AbstractVector{<:Real}, coefa::AbstractVector{<:Real})
    DSP.filt!(data, coefb, coefa, data)
    reverse!(data, dims=1)
    DSP.filt!(data, coefb, coefa, data)
    return reverse!(data, dims=1)
end
# ============================================================================ #
@inline interpolate_channels!(data::AbstractMatrix{<:Real}, ::NullInterpolator) = data

interpolate_channels!(data::AbstractMatrix{<:Real}, ri::RealInterpolator) = interpolate_channels!(data, ri.bad, ri.n)

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
