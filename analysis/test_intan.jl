using Intan, CSV, TensorOps

idir = "/home/scottie/data/intan"

datafile = joinpath(idir, "lowpass.dat")
evtfile = joinpath(idir, "digitalin.dat")

ds = Intan.IntanData(datafile)

ratio = 1 // 30  # resampling ratio

new_fs = Int(ds.fs * ratio)

npre = floor(Int, 0.02 * ds.fs)
npost = floor(Int, 0.3 * ds.fs)

order = channel_order(ds)

# ---------------------------------------------------------------------------- #
# interpolate channel w/ impedance > 1 MOhm
d = CSV.parse(joinpath(idir, "probe_174_first_penetration.csv"), ',', [String,String,String,Bool,Float64,Int,Float64,Float64]);
bad = findall(>(1e6), d["Impedance Magnitude at 1000 Hz (ohms)"][order])
# ---------------------------------------------------------------------------- #
# only pins 1 & 2 are used
seq, onset = Intan.event_sequence(Intan.read_digital_file(evtfile), 1:2)

out = Intan.find_sequence(seq, [2,1,1,1,1,1,1,2])

# 2nd, 4th, and 6th event onset are image onsets, as the sequence is:
# trial_onset(2) -> image_onset (1) -> grey screen (1) -> image_onset (1) -> ...
konset = vcat(map(x->x[[2,4,6]], out)...)
# ---------------------------------------------------------------------------- #
# resample and interpolate bad channels
proc = preprocessor(ratio, bad, 2)

# data = load_and_process(ds, onset[konset], npre, npost, proc);
#
# nbase = floor(Int, 0.02 * new_fs) - 1
# erp = rm_baseline!(mean_3(data), nbase)
