module TargeneCore


if occursin("Intel", Sys.cpu_info()[1].model)
    using MKL
end

using ArgParse
using DataFrames
using CSV
using BGEN
using SnpArrays
using YAML
using Combinatorics
using TMLE
using Arrow
using Serialization
using StableRNGs
using Random
using Distributions
using Statistics
using JLD2
using CairoMakie
using TMLECLI
using HypothesisTests
using OrderedCollections
using Mmap

###############################################################################
###                               INCLUDES                                  ###
###############################################################################

include("utils.jl")
include("confounders.jl")
include("dataset.jl")
include("outputs.jl")
include("inputs_from_estimands.jl")
include("inputs_from_config.jl")
include("estimation_inputs.jl")
include("sieve_variance.jl")
include("cli.jl")

###############################################################################
###                               EXPORTS                                  ###
###############################################################################

export estimation_inputs
export make_dataset
export make_outputs
export filter_chromosome, merge_beds, adapt_flashpca

export get_outcome, get_treatments, get_outcome_extra_covariates, get_confounders, get_all_confounders
export pvalue_or_nan

end
