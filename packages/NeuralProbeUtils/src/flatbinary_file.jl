# ============================================================================ #
struct FlatBinaryFile{P<:AbstractProbe,T}
    probe::P
    filepath::String
    total_channels::Int
    n_sample::Int
    fs::Float64
end

n_channel(bf::FlatBinaryFile) = bf.total_channels
n_sample(bf::FlatBinaryFile) = bf.n_sample
channel_order(bf::FlatBinaryFile) = channel_order(bf.probe)
reference_channels(bf::FlatBinaryFile) = reference_channels(bf.probe)

function memmap(bf::FlatBinaryFile{P,T}) where {P,T}
    return open(bf.filepath, "r") do io
        Mmap.mmap(io, Matrix{T}, (n_channel(bf), n_sample(bf)), grow=false)
    end
end
# ============================================================================ #
