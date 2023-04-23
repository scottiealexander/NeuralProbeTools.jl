push!(LOAD_PATH, joinpath(@__DIR__, "utils"))
push!(LOAD_PATH, joinpath(@__DIR__, "packages"))
using Pkg
Pkg.activate(@__DIR__)
