
"""
    ATEEstimator(target_cond_expectation_estimator, treatment_cond_likelihood_estimator)

# Scope:

Implements the Targeted Minimum Loss-Based estimator of the Average Treatment Effect.
The Average Treatment Effect is defined as: ATE = E[E[Y|T=1, W=w] - E[Y|T=0, W=w]],
where:

- Y is the target variable (Binary)
- T is the treatment variable (Binary)
- W is a confounder variable

The TMLE procedure relies on plugin estimators. Here, the ATE requires both estimators of t,w → E[Y|T=t, W=w]
and w → p(w). The empirical distribution will be used for w → p(w) all along. An estimator of 
t,w → E[Y|T=t, W=w] is then needed. A first estimator is obtained via the provided 
`target_cond_expectation_estimator` argument. Because those initial estimators are not 
tailored to target the ATE they will be biased. The second stage of TMLE consists in optimizing 
bias and variance of the ATE by fluctuating the t,w → E[Y|T=t, W=w] function.

More information can be found about TMLE in "Causal Inference for Observational and Experimental Data"
by Mark J. van der Laan and Sherri Rose.

# Arguments:

- target_cond_expectation_estimator::MLJ.Supervised : Typically a Stack
- treatment_cond_likelihood_estimator::MLJ.Supervised : Typically a Stack

# Examples:

TODO

"""
mutable struct ATEEstimator <: TMLEstimator 
    target_cond_expectation_estimator::MLJ.Supervised
    treatment_cond_likelihood_estimator::MLJ.Supervised
    fluctuation_family::Distribution
end






###############################################################################
## Fluctuation
###############################################################################

function compute_fluctuation(fitted_fluctuator::GeneralizedLinearModel, 
                             target_expectation_mach::Machine, 
                             treatment_likelihood_mach::Machine, 
                             W, 
                             t_target)
    T = (t=float(t_target),)
    X = merge(T, W)
    offset = compute_offset(target_expectation_mach, X)
    cov = compute_covariate(treatment_likelihood_mach, W, T, t_target)
    return  GLM.predict(fitted_fluctuator, reshape(cov, :, 1); offset=offset)
end

###############################################################################
## Fit
###############################################################################


function MLJ.fit(tmle::ATEEstimator, 
                 verbosity::Int, 
                 t::CategoricalVector{Bool}, 
                 W, 
                 y::Union{CategoricalVector{Bool}, Vector{<:Real}})
    n = nrows(y)

    # Convert to NamedTuples
    W = Tables.columntable(W)
    T = (t=float(t),)
    X = merge(T, W)
    t_target = t
    
    # Initial estimate of E[Y|A, W]
    target_expectation_mach = machine(tmle.target_cond_expectation_estimator, X, y)
    fit!(target_expectation_mach, verbosity=verbosity)

    # Estimate of P(A|W)
    treatment_likelihood_mach = machine(tmle.treatment_cond_likelihood_estimator, W, t_target)
    fit!(treatment_likelihood_mach, verbosity=verbosity)
    
    # Fluctuate E[Y|A, W] 
    # on the covariate and the offset 
    offset = compute_offset(target_expectation_mach, X)
    covariate = compute_covariate(treatment_likelihood_mach, W, T, t_target)
    fluctuator = glm(reshape(covariate, :, 1), y, tmle.fluctuation_family; offset=offset)

    # Compute the final estimate tmleATE = 1/n ∑ Fluctuator(t=1, W=w) - Fluctuator(t=0, W=w)
    fluct = (compute_fluctuation(fluctuator, 
                                target_expectation_mach, 
                                treatment_likelihood_mach, 
                                W, 
                                categorical(ones(Bool, n), levels=levels(t_target)))
            - compute_fluctuation(fluctuator, 
                                target_expectation_mach, 
                                treatment_likelihood_mach,
                                W, 
                                categorical(zeros(Bool, n), levels=levels(t_target))))

    estimate = mean(fluct)
    
    # Standard error from the influence curve
    observed_fluct = GLM.predict(fluctuator, reshape(covariate, n, 1); offset=offset)
    inf_curve = covariate .* (float(y) .- observed_fluct) .+ fluct .- estimate

    fitresult = (
        estimate=estimate,
        stderror=sqrt(var(inf_curve)/n),
        target_expectation_mach=target_expectation_mach,
        treatment_likelihood_mach=treatment_likelihood_mach,
        fluctuation=fluctuator
        )

    cache = nothing
    report = nothing
    return fitresult, cache, report

end
