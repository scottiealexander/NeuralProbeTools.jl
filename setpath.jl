push!(LOAD_PATH, joinpath(@__DIR__, "utils"))
push!(LOAD_PATH, joinpath(@__DIR__, "packages"))
push!(LOAD_PATH, joinpath(@__DIR__, "analysis"))
using Pkg
Pkg.activate(@__DIR__)
