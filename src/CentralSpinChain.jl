module CentralSpinChain
__precompile__()

using ITensors
using HDF5
using Plots
using ProgressMeter
using LaTeXStrings

include("model.jl")

export run_simulation

end  # module CentralSpinChain
