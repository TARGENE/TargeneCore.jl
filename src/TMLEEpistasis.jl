module TMLEEpistasis

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

import MLJBase: fit, transform, target_scitype, input_scitype


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
