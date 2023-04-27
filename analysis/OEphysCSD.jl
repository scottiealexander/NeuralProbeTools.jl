module OEphysCSD

using Statistics
using OEphys, TensorOps, CSDView
using NeuralProbeUtils

# basedir = "/home/scottie/data/oephys/recording3"
# apparently bad channels
const REC3_BAD = [14,26,40,42,44,54,56,58,60,70,88,90,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128]
# ============================================================================ #
function get_erp(basedir::AbstractString, bad_channels::AbstractVector{<:Integer}=Int[])

    meta = OEphys.metadata(basedir)

    filepath, nchan, nsample, fs, cont_idx = OEphys.continuous_info(basedir, meta)
    ds = OEphys.OEData(filepath, nchan, nsample, fs)

    # ------------------------------------------------------------------------ #
    # # event times (in seconds) in OE time base
    # evt = OEphys.csd_events(basedir)
    #
    # # sample = round(Int, evt * fs) .+ 1 # sample # at which event occured (index of first sample is usually > 1)
    # # index = sample .- (cont_idx[1] - 1) # convert sample # to index within data file array
    # # or equivalently, subtract 2 from cont_idx[1]...
    # evt_idx = round.(Int, evt * fs) .- (cont_idx[1] - 2)
    # ------------------------------------------------------------------------ #
    onset, dur, lab = OEphys.ttl_events(basedir)
    seqs = OEphys.find_sequence(lab, [2,1,1,1,1,1,1,2])

    # the first 5 1's represent transitions x->white|black, the 6th is white->gray
    kevt = vcat(map(x->x[2:6],seqs)...)
    evt_idx = onset[kevt] .- (cont_idx[1] - 1)
    # ------------------------------------------------------------------------ #

    ratio = 1 // 30  # resampling ratio
    lowcutoff = 0.5  # highpass filter cutoff
    highcutoff = 0.0 # lowpass filter cutoff (0 = omit) (resampling includes lowpass filtering)

    pre = -0.05
    post = 0.25

    # could be Float64(), but this acts as a check that <new_fs> is an integer factor of <fs>
    new_fs = Int(fs * ratio)

    npre = floor(Int, abs(pre) * ds.fs)
    npost = floor(Int, post * ds.fs)

    if isempty(bad_channels)
        # preprocessing does not include interpolation
        proc = preprocessor(ratio, new_fs, lowcutoff, highcutoff)
        data = load_and_process(ds, evt_idx, npre, npost, proc)

        # indices of bad channels from mean erp
        bad_channels = find_bad_channels(mean_3(data), 0.5)

        # inplace interpolation of bad channels
        interpolate_bad_channels!(data, bad_channels, 2)
    else
        proc = preprocessor(ratio, new_fs, lowcutoff, highcutoff, bad_channels, 2)
        data = load_and_process(ds, evt_idx, npre, npost, proc)
        # n = @allocated((data = load_and_process(ds, evt_idx, npre, npost, proc)))
        # total = (post - pre) * fs * size(data, 2) * size(data, 3) * sizeof(Int16)
        # @info("Allocation $(n/2^20), $(total/2^20) $(n/total)")
    end


    nbase = floor(Int, 0.05 * new_fs) - 1
    return rm_baseline!(mean_3(data), nbase), bad_channels
end
# ============================================================================ #
function view_csd(erp, pre, post, title)
    # after OEphys.channel_order is applied to data as it's read in, channels are
    # ordered superficial -> deep (or far-from-probe-tip to near-probe-tip)
    # so no need to reverse (unlike w/ Neuropixel data)
    return CSDView.view(
        erp, (8,3), spacing=OEphys.CONTACT_SPACING,
        xlim=[pre, post], ylim=[3175,0], title=title, rev=false
    )
end
# ============================================================================ #
end
