module GenesInteraction

using DataFrames
using CSV
using SnpArrays
using MLJ
using Distributions
using TMLE
using TOML
using BGEN

###############################################################################
# INCLUDES

include("ukbb.jl")
include("models.jl")
include("stackbuilding.jl")


###############################################################################
# EXPORTS

export UKBBtmleepistasis


end
