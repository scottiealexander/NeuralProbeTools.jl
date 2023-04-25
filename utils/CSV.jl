module CSV

const TypeVec = Vector{DataType}
# ============================================================================ #
mutable struct Row{T<:AbstractString}
    data::T
    delim::Char
end
# ---------------------------------------------------------------------------- #
Base.iterate(r::Row) = iterate(r, 1)
# ---------------------------------------------------------------------------- #
function Base.iterate(r::Row, klast::Int)
    klast < 1 && return nothing

    while isspace(r.data[klast])
        klast += 1
    end

    if r.data[klast] == '"'
        # find the next '"' that is not escaped (preceded by a '\')
        # k = search(r.data, '"', klast+1)
        k = findnext(isequal('"'), r.data, klast+1)
        while r.data[k-1] == '\\'
            # k = search(r.data, '"', k+1)
            k = findnext(isequal('"'), r.data, k+1)
        end

        # k points to the char after the closing '"'
        k += 1

        # knext needs to point to the char after the closing r.delim or 0
        # knext = search(r.data, r.delim, k)
        knext = findnext(isequal(r.delim), r.data, k)
        knext = knext == nothing ? 0 : knext + 1
    else
        # element is not quoted, simply find the next delim
        # k = search(r.data, r.delim, klast)
        k = findnext(isequal(r.delim), r.data, klast)
        if k == nothing
            knext = 0
            k = length(r.data) + 1
        elseif k == length(r.data)
            # string ends with delim (i.e. last element of this row is empty)
            knext = 0
        else
            knext = k + 1
        end
    end
    return SubString(r.data, klast, k-1), knext
end
# ---------------------------------------------------------------------------- #
# Base.done(r::Row, k::Int) = k < 1
Base.IteratorSize(r::Row) = SizeUnknown()
Base.IteratorEltype(r::Row) = HasEltype()
# ---------------------------------------------------------------------------- #
Base.eltype(r::Row{T}) where T<:AbstractString = T
# ---------------------------------------------------------------------------- #
function Base.collect(r::Row)
    out = Vector{String}()
    for x in r
        push!(out, x)
    end
    return out
end
# ============================================================================ #
function parse(ifile::String, delim::Char=',', types::TypeVec=TypeVec())
    return try
        if isfile(ifile)
            return parse_impl(read(ifile, String), delim, types)
        elseif isdir(dirname(ifile))
            error("Input appears to be an invalid file path...")
        else
            return parse_impl(ifile, delim, types)
        end
    catch err
        if isa(err, Base.IOError)
            return parse_impl(ifile, delim, types)
        else
            rethrow(err)
        end
    end
end
# ---------------------------------------------------------------------------- #
function parse_impl(csvstr::String, delim::Char=',', types::TypeVec=TypeVec())
    lines = split(csvstr, r"\n|\r\n|\r", limit=0, keepempty=false)
    hdr = split_line(lines[1], delim)
    if isempty(types)
        types = fill(Any, length(hdr))
    elseif length(types) < length(hdr)
        append!(types, fill(Any, length(hdr) - length(types)))
    end

    data = Dict{String, Vector}(
        (hdr[k] => Vector{types[k]}(undef, length(lines)-1) for k in 1:length(hdr))
    )

    for k = 2:length(lines)
        j = 1
        for item in Row(lines[k], delim)
            data[hdr[j]][k-1] = parse_item(types[j], item)
            j += 1
        end
    end
    return data
end
# ============================================================================ #
isquoted(str::AbstractString) = occursin(r"^\"((?:[^\"\\]|\\.)*)\"$", strip(str))
# ============================================================================ #
dequote(str::AbstractString) = replace(str, r"(?:^\s*\")|(?:\"\s*$)"=>"")
# ============================================================================ #
function split_line(str::AbstractString, delim::Char)
    out = Vector{String}()
    for x in Row(str, delim)
        push!(out, dequote(x))
    end
    return out
end
# ============================================================================ #
function parse_item(::Type{T}, item::AbstractString) where {T<:AbstractFloat}
    return isempty(item) ? NaN : Base.parse(T, strip(item))
end
parse_item(::Type{T}, item::AbstractString) where {T<:Integer} = Base.parse(T, strip(item))
# ============================================================================ #
function parse_item(::Type, item::AbstractString)
    return isquoted(item) ? dequote(item) : strip(item)
end
# ============================================================================ #
end
