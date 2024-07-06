module CentralSpinChain

export run_simulation, plot_results

# 引入相关包
using ITensors
using HDF5
using Plots
using ProgressMeter
using LaTeXStrings

# 引入其他文件
include("model.jl")
include("datagenerate.jl")

end  # module CentralSpinChain
