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
include("random_variants_test.jl")

###############################################################################
###                               EXPORTS                                  ###
###############################################################################

export permutation_tests_tl_inputs
export generate_random_variants_estimands
export tl_inputs_from_actors
export tl_inputs_from_param_files
export filter_chromosome, merge_beds, adapt_flashpca
export tl_inputs

end
