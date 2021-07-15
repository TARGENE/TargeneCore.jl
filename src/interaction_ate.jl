
"""


"""
mutable struct InteractionATEEstimator <: TMLEstimator 
    target_cond_expectation_estimator::MLJ.Supervised
    treatment_cond_likelihood_estimator::MLJ.Supervised
    fluctuation_family::Distribution
end


function MLJ.fit(tmle::InteractionATEEstimator, 
                 verbosity::Int, 
                 T,
                 W, 
                 y::Union{CategoricalVector{Bool}, Vector{<:Real}})
    n = nrows(y)

    # Fit Encoding of the treatment variable
    hot_mach = machine(OneHotEncoder(), T)
    fit!(hot_mach, verbosity=verbosity)

    # Input checks and reformating
    X = combinetotable(hot_mach, t, W)
    
    # Initial estimate of E[Y|A, W]
    target_expectation_mach = machine(tmle.target_cond_expectation_estimator, X, y)
    fit!(target_expectation_mach, verbosity=verbosity)

    # Estimate of P(A|W)
    treatment_likelihood_mach = machine(tmle.treatment_cond_likelihood_estimator, W, t)
    fit!(treatment_likelihood_mach, verbosity=verbosity)
end