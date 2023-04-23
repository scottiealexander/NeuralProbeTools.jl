module CSD

using TensorOps, LinearAlgebra

@enum CSDMethod D2 iCSD ERP
# ============================================================================ #
# NOTE: erp is assumed to be channels x time
function csd(erp::AbstractMatrix{<:Real}, smoothing::Tuple{Int,Int}=(8,0), spacing::Real=1, method::CSDMethod=D2)

    if sum(smoothing) > 0
        # using KernelFactors.gaussian() seems to introduce some kind of ringing
        # that is not noticable in the smoothed ERP but dominates the CSD
        # IIRGaussian() does not have this problem, no idea why...
        d = smooth(erp, smoothing) # uses IIRGaussian by default
    else
        d = erp
    end

    nrow = size(d, 1)

    if method == D2
        # 2nd derivative matrix operator (for a single contact column)
        ddp = -d2(nrow, spacing)
    elseif method == iCSD
        # Peterson et al. 2006 J Neurosci Meth inverse CSD method
        ddp = icsd(nrow, spacing)
    else
        ddp = Matrix(I, nrow, nrow)
    end

    return ddp * d
end
# ============================================================================ #
function icsd(n::Integer, h::Real)
    conductivity = 0.3
    radius = 5 * h
    z = h:h:h*n # electrode positions
    F = zeros(n, n)
    for j = 1:n
        for k = 1:n
            F[j,k] = (h/(2*conductivity)) * (sqrt((z[j]-z[k])^2+(radius^2))-abs(z[j]-z[k]));
        end
    end
    return pinv(F)
end
# ============================================================================ #
"2nd derivative matrix operator"
function d2(n::Integer, h::Real)
    out = zeros(n-2, n)
    for k = 1:(n-2)
        for j = 1:n
            if k == j-1
                out[k,j] = -2/h^2
            elseif abs(k-j+1) == 1
                out[k,j] = 1/h^2
            end
        end
    end
    return out
end
# ============================================================================ #
end
