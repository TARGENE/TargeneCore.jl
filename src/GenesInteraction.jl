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

export ATEEstimator, fit

# #############################################################################
# INCLUDES
# #############################################################################

include("tmle.jl")


end
