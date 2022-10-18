module TargeneCore

if occursin("Intel", Sys.cpu_info()[1].model)
    using MKL
end
using DataFrames
using CSV
using TMLE
using BGEN
using JLD2
using SnpArrays
using Mmap
using HypothesisTests
using YAML
using Combinatorics
using TargetedEstimation

###############################################################################
# INCLUDES

include("confounders.jl")
include("sieve_plateau.jl")
include(joinpath("tmle_inputs", "tmle_inputs.jl"))
include(joinpath("tmle_inputs", "from_actors.jl"))
include(joinpath("tmle_inputs", "from_param_files.jl"))


###############################################################################
# EXPORTS

export filter_ukb_genetic_file, merge_beds, adapt_flashpca
export sieve_variance_plateau
export tmle_inputs

end
