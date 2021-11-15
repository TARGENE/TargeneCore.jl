module GenesInteraction

using DataFrames
using CSV
using SnpArrays
using MLJ
using Distributions
using TMLE
using TOML
using BGEN
using HighlyAdaptiveLasso

###############################################################################
# GENERIC FUNCTIONS

is_binary(y) = sort(unique(y)) == [0, 1]


###############################################################################
# INCLUDES

include("ukbb.jl")
include("models.jl")
include("stackbuilding.jl")


###############################################################################
# EXPORTS

export TMLEEpistasisUKBB


end
