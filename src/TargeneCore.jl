module TargeneCore

using DataFrames
using CSV
using BGEN
using SnpArrays
using Mmap
using YAML
using Combinatorics
using TMLE
using Arrow
using Serialization
using StableRNGs
using Random
using JSON
using HTTP
using Statistics
using JLD2

###############################################################################
###                               INCLUDES                                  ###
###############################################################################

include("utils.jl")
include("confounders.jl")
include(joinpath("tl_inputs", "tl_inputs.jl"))
include(joinpath("tl_inputs", "from_actors.jl"))
include(joinpath("tl_inputs", "from_param_files.jl"))
include(joinpath("tl_inputs", "allele_independent_estimands.jl"))
include(joinpath("tl_inputs", "permutation_test.jl"))
include(joinpath("tl_inputs", "random_variants_test.jl"))

###############################################################################
###                               EXPORTS                                  ###
###############################################################################

export generate_permutation_parameters_and_dataset, generate_random_variants_parameters_and_dataset
export filter_chromosome, merge_beds, adapt_flashpca
export tl_inputs

end
