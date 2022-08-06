module TargeneCore

using Base.Sys

if occursin("Intel", cpu_info()[1].model)
    using MKL
end
using DataFrames
using MLJBase
using CSV
using TMLE
using TOML
using BGEN
using Serialization
using JLD2
using SnpArrays
using Mmap
using HypothesisTests
using DelimitedFiles
using YAML
using Combinatorics

import MLJBase: fit, transform, target_scitype, input_scitype

###############################################################################
# INCLUDES

include("utils.jl")
include("confounders.jl")
include("sieve_plateau.jl")
include("summary.jl")
include("tmle_inputs.jl")

###############################################################################
# EXPORTS

export filter_chromosome, merge_beds, adapt_flashpca
export sieve_variance_plateau
export build_summary
export tmle_inputs

end
