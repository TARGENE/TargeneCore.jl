module GenesInteraction

using Distributions
using CategoricalArrays
using MLJ
using GLM: glm, GeneralizedLinearModel
using GLM: predict as predict_glm

import MLJ.fit!

# #############################################################################
# EXPORTS
# #############################################################################

export ATEEstimator, fit!

# #############################################################################
# INCLUDES
# #############################################################################

include("tmle.jl")


end
