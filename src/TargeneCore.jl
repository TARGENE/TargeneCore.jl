module TargeneCore

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
using TargetedEstimation

###############################################################################
###                               INCLUDES                                  ###
###############################################################################

include("utils.jl")
include("confounders.jl")
include("dataset.jl")
include("plots.jl")
include(joinpath("tl_inputs", "tl_inputs.jl"))
include(joinpath("tl_inputs", "from_param_files.jl"))
include(joinpath("tl_inputs", "allele_independent_estimands.jl"))
include("estimation_inputs.jl")
include("cli.jl")

###############################################################################
###                               EXPORTS                                  ###
###############################################################################

export tl_inputs
export tl_inputs_from_param_files

export generate_dataset
export generate_summary_plots
export filter_chromosome, merge_beds, adapt_flashpca

end
