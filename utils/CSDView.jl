module CSDView

using CSD, InteractiveImage, PyPlot, TensorOps

# ============================================================================ #
function view(erp::AbstractMatrix{<:Real}, smoothing::Tuple{Int,Int};
    method::CSD.CSDMethod=CSD.D2, spacing::Real=1.0,
    ylim::Vector{<:Real}=[size(erp,2),1], xlim::Vector{<:Real}=[1,size(erp,1)],
    rev::Bool=true, erp_smoothing::Tuple{Int,Int}=(0,2), title::AbstractString="")

    csd = CSD.csd(erp', smoothing, spacing, method)

    h, ax = subplots(1,2)
    h.set_size_inches((12,6))

    show_image(
        (sum(erp_smoothing) > 0 ? smooth(erp, erp_smoothing) : erp)',
        xlim, ylim, "ERP", title, rev, ax[1]
    )
    show_image(csd, xlim, ylim, "CSD", title, rev, ax[2])

    return h, ax
end
# ============================================================================ #
function show_image(im, xlim, ylim, type, title, rev, ax=nothing)

    him, ax = InteractiveImage.display(
        rev ? reverse(im, dims=1) : im,
        xlim, ylim, ax=ax
    )

    ax.plot([0,0], ax.get_ylim(), "--", linewidth=2, color="white")

    ax.set_xlabel("Time (seconds, 0 = stim onset)", fontsize=14)
    ax.set_ylabel("Relative depth (\$\\mu m\$)", fontsize=14)

    ax.set_title("$(type): $(title)", fontsize=18)
    return him, ax
end
# ============================================================================ #
end
