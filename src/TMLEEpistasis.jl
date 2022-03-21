module TMLEEpistasis

using MKL
using DataFrames
using MLJBase
using CSV
using TMLE
using TOML
using BGEN
using HighlyAdaptiveLasso
using EvoTrees
using MLJModels
using MLJLinearModels
using Serialization
using JLD2
using SnpArrays
using Mmap
using HypothesisTests

import MLJBase: fit, transform, target_scitype, input_scitype


###############################################################################
# INCLUDES

include("ukbb.jl")
include("models.jl")
include("estimators.jl")
include("utils.jl")
include("asb_snps.jl")
include("queries_generation.jl")
include("confounders.jl")
include("phenotypes.jl")
include("grm.jl")
include("sieve_plateau.jl")
include("summary.jl")

###############################################################################
# EXPORTS

export UKBBVariantRun
export generate_queries
export filter_asb
export filter_chromosome, merge_beds, adapt_flashpca
export prepare_phenotypes, tmle_phenotypes_batches
export sieve_variance_plateau
export build_summary

end
