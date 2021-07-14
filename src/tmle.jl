

abstract type TMLEstimator end

pvalue(tmle::TMLEstimator) = 2*(1 - cdf(Normal(0, 1), abs(tmle.estimate/tmle.std_error)))

confint(tmle::TMLEstimator) = (tmle.estimate - 1.96tmle.stderror, tmle.estimate + 1.96tmle.stderror)


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
## Helper functions
###############################################################################


function combinetotable(hot_mach::Machine, t, W)
    thot = MLJ.transform(hot_mach, (t=t,))
    return merge(thot, columntable(W))
end


logit(X) = log.(X ./ (1 .- X))

"""
Hack into GLM to compute deviance on y a real
"""
function GLM.devresid(::Bernoulli, y::Vector{<:Real}, μ::Real)
    return -2*(y*log(μ) + (1-y)*log1p(-μ))
end

"""
Remove default check for y to be binary
"""
GLM.checky(y, d::Bernoulli) = nothing

###############################################################################
## Offset and covariate
###############################################################################

function compute_offset(fitted_mach::Machine{<:Probabilistic}, X)
    # The machine is an estimate of a probability distribution
    # In the binary case, the expectation is assumed to be the probability of the second class
    expectation = MLJ.predict(fitted_mach, X).prob_given_ref[2]
    return logit(expectation)
end


function compute_offset(fitted_mach::Machine{<:Deterministic}, X)
    return MLJ.predict(fitted_mach, X)
end


function compute_covariate(fitted_mach::Machine, W, t)
    # tpred is an estimate of a probability distribution
    # we need to extract the observed likelihood out of it
    tint = convert(Vector{Int}, t)
    tpred = MLJ.predict(fitted_mach, W)
    likelihood_fn = pdf(tpred, levels(first(tpred)))
    likelihood = [x[tint[i] + 1] for (i, x) in enumerate(eachrow(likelihood_fn))]
    # truncate predictions, is this really necessary/suitable?
    likelihood = min.(0.995, max.(0.005, likelihood))
    return (2tint .- 1) ./ likelihood
end

###############################################################################
## Fluctuation
###############################################################################

function compute_fluctuation(fitted_fluctuator::GeneralizedLinearModel, 
                             target_expectation_mach::Machine, 
                             treatment_likelihood_mach::Machine, 
                             hot_mach::Machine,
                             W, 
                             t::CategoricalVector{Bool})
    X = combinetotable(hot_mach, t, W)
    offset = compute_offset(target_expectation_mach, X)
    cov = compute_covariate(treatment_likelihood_mach, W, t)
    return  GLM.predict(fitted_fluctuator, reshape(cov, :, 1); offset=offset)
end

###############################################################################
## Fit
###############################################################################


function MLJ.fit(tmle::TMLEstimator, 
                 verbosity::Int, 
                 t::CategoricalVector{Bool}, 
                 W, 
                 y::Union{CategoricalVector{Bool}, Vector{<:Real}})
    n = nrows(y)

    # Fit Encoding of the treatment variable
    hot_mach = machine(OneHotEncoder(), (t=t,))
    fit!(hot_mach, verbosity=verbosity)

    # Input checks and reformating
    X = combinetotable(hot_mach, t, W)
    
    # Initial estimate of E[Y|A, W]
    target_expectation_mach = machine(tmle.target_cond_expectation_estimator, X, y)
    fit!(target_expectation_mach, verbosity=verbosity)

    # Estimate of P(A|W)
    treatment_likelihood_mach = machine(tmle.treatment_cond_likelihood_estimator, W, t)
    fit!(treatment_likelihood_mach, verbosity=verbosity)
    
    # Fluctuate E[Y|A, W] 
    # on the covariate and the offset 
    offset = compute_offset(target_expectation_mach, X)
    covariate = compute_covariate(treatment_likelihood_mach, W, t)
    fluctuator = glm(reshape(covariate, :, 1), y, tmle.fluctuation_family; offset=offset)

    # Compute the final estimate tmleATE = 1/n ∑ Fluctuator(t=1, W=w) - Fluctuator(t=0, W=w)
    fluct_treatment_true = compute_fluctuation(fluctuator, 
                                               target_expectation_mach, 
                                               treatment_likelihood_mach, 
                                               hot_mach,
                                               W, 
                                               categorical(ones(Bool, n), levels=[false, true]))
    fluct_treatment_false = compute_fluctuation(fluctuator, 
                                                target_expectation_mach, 
                                                treatment_likelihood_mach,
                                                hot_mach,
                                                W, 
                                                categorical(zeros(Bool, n), levels=[false, true]))

    estimate = mean(fluct_treatment_true .- fluct_treatment_false)
    
    # Standard error from the influence curve
    observed_fluct = GLM.predict(fluctuator, reshape(covariate, n, 1); offset=offset)
    inf_curve = covariate .* (float(y) .- observed_fluct) .+ fluct_treatment_true .- fluct_treatment_false .- estimate

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
