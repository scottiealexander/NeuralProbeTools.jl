# ============================================================================ #
# the "probe" interface, defined as sub-types of AbstractProbe
abstract type AbstractProbe end

# **MUST** define the following methods on an instance
"returns the number of channels on the physical probe"
n_channel(::AbstractProbe) = error("Not implemented")

"returns a permutation Vector{Int} that puts the channels in a superficial->deep ordering"
channel_order(::AbstractProbe) = error("Not implemented")

"return the vertical channel pitch in meters"
vertical_pitch(::AbstractProbe) = error("Not implemented")

"returns a Vector{Int} of reference channel indices"
reference_channels(::AbstractProbe) = Int[]

"""
`pos = channel_positions(probe::AbstractProbe)`

returns an Nx2 Float64 matrix where each row is the (x,y) position of the
corresponding channel in the binary file, i.e.:

`pos[1,1]` -> the horizontal position of the first channel listed in the binary file

`pos[1,2]` -> the vertical position of the first channel listed in the binary file

**NOTE**: for vertical positions, larger numbers mean closer to the probe tip
which **usually** means deeper in the brain
"""
channel_positions(::AbstractProbe) = error("Not implemented")

# ============================================================================ #
struct Neuropixel3A <: AbstractProbe end
n_channel(::Neuropixel3A) = 384
channel_order(::Neuropixel3A) = 384:-1:1 # flip channel order to be shallow -> deep, omit sync channel # 385
vertical_pitch(::Neuropixel3A) = 0.000020
reference_channels(::Neuropixel3A) = [37, 76, 113, 152, 189, 228, 265, 304, 341, 380]
channel_positions(::Neuropixel3A) = read_channel_positions(joinpath(@__DIR__, "..", "data", "neuropixel3a.f64.dat"))
# ============================================================================ #
struct DBCDeepArray <: AbstractProbe end
n_channel(::DBCDeepArray) = 128
channel_order(::DBCDeepArray) = reshape(vcat((1:64)',(128:-1:65)'),128) # put channels into shallow (1) -> deep (65) ordering
vertical_pitch(::DBCDeepArray) = 0.000025
channel_positions(::DBCDeepArray) = read_channel_positions(joinpath(@__DIR__, "..", "data", "dbcdeeparray.f64.dat"))
# ============================================================================ #
# NOTE: channels position files hold the x and y corrdinates of each channel
# in an N x 2 matrix (call it `pos`) such that `pos[1,1]` and `pos[1,2]` are the
# x and y corrdinates (respectivly) of the first channel listed in the data file
# NOTE: that depths (y-coords) are reltive such that larger number mean "closer
# to the probe tip" and thus usually deeper in the brain
function read_channel_positions(ifile::AbstractString)
    return open(ifile, "r") do io
        dims = DataChunks.read_dims(io)
        pos = zeros(Float64, dims[1], dims[2])
        read!(io, pos)
        return pos
    end
end
# ============================================================================ #
