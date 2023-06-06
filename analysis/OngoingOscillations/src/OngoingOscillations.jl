module OngoingOscillations

using DSP
using NeuralProbeUtils, TensorOps, GLX

import NeuralProbeUtils, Base

const DATADIR = "E:\\SGL_DATA\\"

include("./preprocessing.jl")
# ============================================================================ #

# ============================================================================ #
function get_npxdir(run::AbstractString; artefact_corrected::Bool=false, datadir::AbstractString=DATADIR)
    suffix = artefact_corrected ? "artefact_correct" : ""
    return get_npxdir(run, suffix, datadir)
end
# ---------------------------------------------------------------------------- #
function get_npxdir(runname::AbstractString, suffix::AbstractString, datadir::AbstractString=DATADIR)

    names = filter!(readdir(datadir)) do name
        occursin(runname, name)
    end

    if isempty(names)
        error("Failed to locate \"$(runname)\" in \"$(datadir)\"")
    elseif length(names) > 1
        error("Given runname \"$(runname)\" is NOT unique in \"$(datadir)\"")
    end

    return joinpath(datadir, names[1], suffix)
end
# ============================================================================ #
function run(rundir::AbstractString)

    ratio = 2//5
    lowcutoff = 0.0
    highcutoff = 120.0
    pre = -0.05
    post = 1.0
    window_dur = 1.0

    file = GLX.spikeglx_data(rundir, false)

    new_fs = Int(ratio * file.fs)

    # pre-trial tve
    window_pre = round(Int, (abs(pre) + window_dur) * fs)
    window_post = round(Int, pre * fs)

    window_proc = NpxProcessor(ratio, new_fs, lowcutoff, highcutoff, window_pre, window_post)



    # npre = round(Int, abs(pre) * file.fs)
    # npost = round(Int, post * file.fs)
    #
    # trial_proc = NpxProcessor(ratio, new_fs, lowcutoff, highcutoff, npre, npost)

    # # load event information
    #
    # # load erp data
    # data = load_and_process(file, onset, npre, npost, trial_proc)

end
# ============================================================================ #
function pretrial_tve(file::NeuralProbeUtils.FlatBinaryFile, p::NpxProcessor,
    window_dur::Real, pre::Real, post::Real)


    proc = NpxProcessor(p.ratio, )

end
# ============================================================================ #
end
