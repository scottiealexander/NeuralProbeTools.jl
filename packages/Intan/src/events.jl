function event_sequence(evt::Vector{IOType}, bits::AbstractVector=1:16)
    seq = Vector{Int}(undef, 0)
    onset = Vector{Int}(undef, 0)

    for k in bits
        idx, _ = event_indices(evt, k)
        append!(seq, fill(k, length(idx)))
        append!(onset, idx)
    end

    ks = sortperm(onset)

    return seq[ks], onset[ks]
end

function find_sequence(seq::Vector{Int}, target::AbstractVector{<:Integer})
    out = Vector{UnitRange{Int}}(undef, 0)
    len = length(target)
    k = 1
    while k <= (length(seq) - len) + 1
        idx = k:k+len-1
        if seq[idx] == target
            push!(out, idx)#onset[idx])
            k += len
        else
            k += 1
        end
    end
    return out
end

#=

using Intan, GLXUtils, TensorOps, Statistics, CSV, PyPlot

intandir = "/home/scottie/data/intan"
datafile = joinpath(intandir, "lowpass.dat");
evtfile = joinpath(intandir, "digitalin.dat");

fs = Intan.read_sample_rate(datafile)

tmp = Intan.read_digital_file(evtfile);

# t = Intan.read_time_file(joinpath(intandir, "time.dat"));

seq, onset = Intan.event_sequence(tmp, 1:2);

out = Intan.find_sequence(seq, [2,1,1,1,1,1,1,2]);

# 2nd, 4th, and 6th event onset are image onsets, as the sequence is:
# trial_onset(2) -> image_onset (1) -> grey screen (1) -> image_onset (1) -> ...
konset = vcat(map(x->x[[2,4,6]], out)...)

# NOTE: when plotting using <t> as the x variable, using Float64 event times
# appears to produce an incorrect results (when plotted)
evt = (onset[konset] .- 1) ./ fs # <- keep times in Float32 (appears correct) as <fs> is Float32
# evt = t[onset[konset]]

pre, post = -0.02, 0.18
@time data, _ = erp_serial(Intan.erp_access(datafile), evt[1:100], pre, post);

d = CSV.parse(joinpath(intandir, "probe_174_first_penetration.csv"), ',', [String,String,String,Bool,Float64,Int,Float64,Float64]);

bad = findall(>(1e6), d["Impedance Magnitude at 1000 Hz (ohms)"])

@time Intan.interpolate_channels!(data, bad, 2);

mn = rm_baseline!(mean_3(data), 600);

figure(); imshow(mn', extent=[pre, post, 1, 128]); gca().set_aspect(abs(post - pre) / 127)

figure(); plot(1:128, apply_n(std, mn, 1))

=#
