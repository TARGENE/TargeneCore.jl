module GenesInteraction

using Distributions: expectation
using Tables
using Distributions
using CategoricalArrays
using MLJ
using GLM

import MLJ.fit
import GLM.checky
import GLM.devresid

# #############################################################################
# EXPORTS
# #############################################################################

export ATEEstimator, InteractionATEEstimator
export fit

# #############################################################################
# INCLUDES
# #############################################################################

include("utils.jl")
include("tmle.jl")
include("ate.jl")
include("interaction_ate.jl")

end
