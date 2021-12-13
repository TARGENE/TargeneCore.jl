module TMLEEpistasis

using DataFrames
using CSV
using SnpArrays
using MLJ
using Distributions
using TMLE
using TOML
using BGEN
using HighlyAdaptiveLasso

import MLJ: fit, transform, target_scitype

###############################################################################
# GENERIC FUNCTIONS

is_binary(y) = sort(unique(y)) == [0, 1]


###############################################################################
# INCLUDES

include("ukbb.jl")
include("models.jl")
include("estimators.jl")
include("utils.jl")

###############################################################################
# EXPORTS

export UKBBVariantRun

end
