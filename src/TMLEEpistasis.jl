module TMLEEpistasis

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

import MLJBase: fit, transform, target_scitype, input_scitype


###############################################################################
# INCLUDES

include("utils.jl")
include("genotypes_and_queries.jl")
include("confounders.jl")
include("phenotypes.jl")
include("grm.jl")
include("sieve_plateau.jl")
include("summary.jl")
include("tmle_inputs.jl")

###############################################################################
# EXPORTS

export build_genotypes_and_queries
export filter_chromosome, merge_beds, adapt_flashpca
export prepare_phenotypes, tmle_phenotypes_batches
export sieve_variance_plateau
export build_summary
export finalize_tmle_inputs

end
