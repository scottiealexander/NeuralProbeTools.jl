module OEphysCSD

using Statistics
using OEphys, TensorOps, InteractiveImage
using NeuralProbeUtils

# basedir = "/home/scottie/data/oephys/recording3"
# apparently bad channels
const REC3_BAD = [14,26,40,42,44,54,56,58,60,70,88,90,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128]
# ============================================================================ #
"""
`h, erp, bad_channels = OEphysCSD.run(basedir::String; <keyword arguments>)`

Calculate and plot CSD for DeepArray data recorded with OpenEphys

# Inputs:
* `basedir` - the full path to an OpenEphys recording's base directory
                where the \"continuous\" and \"events\" folders can be found
# Options:
* `pre::Real`                 - [-0.05] relative pre-stimulus baseline duration in SECONDS
* `post::Real`                - [0.25] post-stimulus period duration in SECONDS
* `bad_channels::Vector{Int}` - [<auto>] indices of channels to interpolate
* `resample::Rational{Int}`   - [1//30] resampling ratio
* `lowcutoff::Real`           - [0.5] frequency in HZ of highpass filter cutoff(set to 0 to omit)
* `highcutoff::Real`          - [0.0] frequency in HZ of lowpass filter cutoff (set to 0 to omit)
* `csd_smoothing::Tuple{Int,Int}` - [(8,3)] std (space,time) of Gaussian smoothing kernel for CSD in PIXELS
* `erp_smoothing::Tuple{Int,Int}` - [(2,0)] std (space,time) of Gaussian smoothing kernel for ERP in PIXELS
* `title::String`             - [\"\"] the name of the recording (e.g. \"recording 3\")

# Outputs:
* `h` - the handle to the created figure
* `erp` - a `time x channels` ERP matrix
* `bad_channels` - a list of indices of \"bad channels\", when the input option
                     `bad_channels` is not set, bad channels will be automatically
                     identified and the result is output for future calls

# Examples
```
h, erp, bad = OEphysCSD.run(\"/home/user/data/deep-array/recording3\";
                pre=-0.02, post=0.3, resample=1//30, lowcutoff=1.0,
                title="recording 3", erp_smoothing=(2,2)
            )

# reuse bad channels that were automatically identified from recording 3
OEphysCSD.run(\"/home/user/data/deep-array/recording4\", bad_channels=bad, title="recording 4")
```

"""
function run(basedir::AbstractString; pre::Real=-0.05, post::Real=0.25,
    resample::Rational{Int}=1//30, lowcutoff::Real=0.5, highcutoff::Real=0.0,
    bad_channels::AbstractVector{<:Integer}=Int[], csd_smoothing::Tuple{Int,Int}=(8,3),
    args...)
    # ------------------------------------------------------------------------ #
    # # event times (in seconds) in OE time base
    # evt = OEphys.csd_events(basedir)
    #
    # # sample = round(Int, evt * fs) .+ 1 # sample # at which event occured (index of first sample is usually > 1)
    # # index = sample .- (cont_idx[1] - 1) # convert sample # to index within data file array
    # # or equivalently, subtract 2 from cont_idx[1]...
    # evt_idx = round.(Int, evt * fs) .- (cont_idx[1] - 2)
    # ------------------------------------------------------------------------ #
    probe = DBCDeepArray()
    file, onset, dur, lab = OEphys.openephys_data(basedir, probe)

    seqs = OEphys.find_sequence(lab, [2,1,1,1,1,1,1,2])

    # the first 5 1's represent transitions x->white|black, the 6th is white->gray
    kevt = vcat(map(x->x[2:6],seqs)...)
    evt_idx = onset[kevt]

    @info("$(length(evt_idx)) CSD transitions located")
    # ------------------------------------------------------------------------ #

    # could be Float64(), but this acts as a check that <new_fs> is an integer factor of <fs>
    new_fs = Int(file.fs * resample)

    npre = floor(Int, abs(pre) * file.fs)
    npost = floor(Int, post * file.fs)

    if isempty(bad_channels)
        # preprocessing does not include interpolation
        proc = Preprocessor(Resampler(resample), Filterer(new_fs, lowcutoff, highcutoff))
        data = load_and_process(file, evt_idx, npre, npost, proc)

        # indices of bad channels from mean erp
        bad_channels = find_bad_channels(mean_3(data), 0.5)

        # inplace interpolation of bad channels
        interpolate_bad_channels!(data, bad_channels, 2)
    else
        proc = Preprocessor(
            Resampler(resample),
            Filterer(new_fs, lowcutoff, highcutoff),
            Interpolator(bad_channels, 2)
        )
        data = load_and_process(file, evt_idx, npre, npost, proc)
        # n = @allocated((data = load_and_process(file, evt_idx, npre, npost, proc)))
        # total = ((post - pre) * file.fs) * size(data, 2) * size(data, 3) * sizeof(Int16)
        # @info("Allocation $(n/2^20), $(total/2^20) $(n/total)")
    end

    nbase = floor(Int, abs(pre) * new_fs) - 1
    erp = rm_baseline!(mean_3(data), nbase)

    # max and min y-position of channels
    yl = extrema(channel_positions(probe)[:,2])

    h, ax = InteractiveImage.csd_view(erp, csd_smoothing;
        spacing=vertical_pitch(probe), xlim=[pre, post], ylim=[yl[2],yl[1]],
        rev=false, args...)

    return h, erp, bad_channels
end
# ============================================================================ #
function view_csd(erp, pre, post; args...)
    # after OEphys.channel_order is applied to data as it's read in, channels are
    # ordered superficial -> deep (or far-from-probe-tip to near-probe-tip)
    # so no need to reverse (unlike w/ Neuropixel data)
    return InteractiveImage.csd_view(
        erp, (8,3); spacing=OEphys.CONTACT_SPACING,
        xlim=[pre, post], ylim=[3175,0], rev=false, args...
    )
end
# ============================================================================ #
end
