module OngoingOscillations

import NeuralProbeUtils, Base
using NeuralProbeUtils, TensorOps

include("./preprocessing.jl")
# ============================================================================ #
struct TVEProcessor <: AbstractProcessor
    pre::NeuropixelPreprocessor
    freq::Vector{Float64}
    fs::Float64
    ncycle::Int
    sigma::Int
end

Base.copy(p::TVEProcessor) = TVEProcessor(copy(pre), freq, fs, ncycle, sigma)
output_size(p::TVEProcessor, len::Integer, nchan::Integer) = output_size(p.pre, len, nchan)
output_type(::TVEProcessor) = output_type(p.pre)

function process!(out::AbstractMatrix{<:Real}, inp::AbstractMatrix{<:Real}, p::TVEProcessor)
    process!(out, inp, p.pre)
    frequency_envelope!(out, p.freq, p.fs, p.ncycle. p.sigma)
    return out
end
# ============================================================================ #
function run()

end
# ============================================================================ #
end
