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


###############################################################################
# INCLUDES

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
