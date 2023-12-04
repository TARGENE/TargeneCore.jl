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

###############################################################################
###                               INCLUDES                                  ###
###############################################################################

include("confounders.jl")
include(joinpath("tmle_inputs", "tmle_inputs.jl"))
include(joinpath("tmle_inputs", "from_actors.jl"))
include(joinpath("tmle_inputs", "from_param_files.jl"))


###############################################################################
###                               EXPORTS                                  ###
###############################################################################

export filter_chromosome, merge_beds, adapt_flashpca
export tmle_inputs

end
