module InteractiveImage

using PyPlot, PyCall, TensorOps

include("./csd_view.jl")

# ============================================================================ #
function test()
    IMG = randn(200,200)
    xl = 0.001:0.001:0.2
    yl = reverse(1:200)
    display(IMG, xl, yl)
end
# ============================================================================ #
function display(IMG::AbstractMatrix{<:Real}, xl::AbstractVector{<:Real}, yl::AbstractVector{<:Real}; ax=nothing, smoothing::Tuple{Int,Int}=(0,0))

    if ax == nothing
        h = figure()
        ax = PyPlot.axes()
        h.set_size_inches(6, 8, forward=true)
    else
        h = ax.figure
    end

    xmn, xmx = xl[1], xl[end]
    ymn, ymx = yl[1], yl[end]

    if any(>(0), smoothing)
        IMG = smooth(IMG, KernelFactors.IIRGaussian(smoothing))
    end

    him = ax.imshow(IMG, extent=[xmn, xmx, ymn, ymx])

    aspect = abs(xmx - xmn) / abs(ymx - ymn)
    ax.set_aspect(aspect)

    clim = him.get_clim()

    function reset(evt)
        evt.inaxes.images[1].set_clim(clim)
        evt.inaxes.images[1].set_interpolation("none")
    end

    function scroll(evt)
        cl = evt.inaxes.images[1].get_clim()
        evt.inaxes.images[1].set_clim(cl .* (1.0 + (0.02 * evt.step)))
    end

    function key_press(evt)
        if evt.key == "ctrl+w"
            close(h)
        elseif evt.key == "i"
            evt.inaxes.images[1].set_interpolation("bicubic")
        elseif evt.key == "u"
            evt.inaxes.images[1].set_interpolation("none")
        elseif evt.key == "x"
            reset(evt)
        end
    end

    h.canvas.mpl_connect("key_press_event", key_press)
    h.canvas.mpl_connect("scroll_event", scroll)

    return him, ax
end
# ============================================================================ #
end
