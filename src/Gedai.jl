module Gedai

using MAT
using LinearAlgebra
using PosDefManifold
using Statistics
using Random
using DSP
using ThreadPinning
using PrecompileSignatures: @precompile_signatures

# Colors
const font1color = "\x1b[38;5;111m" # metal blue
const font2color = "\x1b[38;5;87m" # cyan
const font3color = "\x1b[38;5;71m" # EEGPlot green
const fontgrey = "\x1b[38;5;249m"
const fontwhite = "\x1b[37m"

# Paths
const leadfield_cov_path = joinpath(@__DIR__, "io", "LF_4_covmat.mat")

# Types
const SymOrHerm{T,S} = Union{Symmetric{T,S}, Hermitian{T,S}} where T<:Real
const Chols = Tuple{LowerTriangular, UpperTriangular}
const Whits = Tuple{SymOrHerm, SymOrHerm}

include("miscellaneous.jl")
include("sensai.jl")
include("broadband.jl")
include("wavelets.jl")
include(joinpath("io", "simple_import.jl"))
include("api.jl")


export  denoise,                # from api.jl
        read_example_data,      # from io/simple_import.jl
        precompute,             # from miscellaneous.jl
        refcov,                 # "
        regularize              # "

@precompile_signatures(Gedai)

end # module

