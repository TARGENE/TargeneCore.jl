

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
    fitted::Bool
    estimate
    stderror
end


function ATEEstimator(target_cond_expectation_estimator::MLJ.Supervised,
                      treatment_cond_likelihood_estimator::MLJ.Supervised)
    return ATEEstimator(target_cond_expectation_estimator,
                        treatment_cond_likelihood_estimator,
                        false,
                        nothing,
                        nothing)
end


function reformat(t::CategoricalVector{Bool}, W, y::CategoricalVector)
    # Ensure target is Binary
    length(levels(y)) == 2 || 
        throw(ArgumentError("y should be a vector of binary values"))
    
    # The treatment will be used both as an input and a target
    tint = convert(Vector{Int}, t)

    X = hcat(tint, MLJ.matrix(W))
    return (X, tint, t, W, y)
end


function compute_offset(fitted_mach::Machine, X)
    # The machine is an estimate of a probability distribution
    # In the binary case, the expectation is assumed to be the probability of the second class
    expectation = predict(fitted_mach, X).prob_given_ref[2]
    return log.(expectation ./ (1 .- expectation))
end


function compute_covariate(fitted_mach::Machine, W, t::Vector{<:Real})
    # tpred is an estimate of a probability distribution
    # we need to extract the observed likelihood out of it
    tpred = predict(fitted_mach, W)
    likelihood_fn = pdf(tpred, levels(first(tpred)))
    likelihood = [x[t[i] + 1] for (i, x) in enumerate(eachrow(likelihood_fn))]
    return (2t .- 1) ./ likelihood
end


function compute_fluctuation(fitted_fluctuator::GeneralizedLinearModel, 
                             target_expectation_mach::Machine, 
                             treatment_likelihood_mach::Machine, 
                             W, 
                             t::Vector{<:Real})
    X = hcat(t, W)
    offset = compute_offset(target_expectation_mach, X)
    cov = compute_covariate(treatment_likelihood_mach, W, t)
    return  predict_glm(fitted_fluctuator, reshape(cov, :, 1); offset=offset)
end


function MLJ.fit!(tmle::ATEEstimator, t::CategoricalVector{Bool}, W, y::CategoricalVector)
    n = nrows(y)
    X, tint, t, W, y = reformat(t, W, y)
    
    # Initial estimate of E[Y|A, W]
    target_expectation_mach = machine(tmle.target_cond_expectation_estimator, X, y)
    fit!(target_expectation_mach)

    # Estimate of P(A|W)
    treatment_likelihood_mach = machine(tmle.treatment_cond_likelihood_estimator, W, t)
    fit!(treatment_likelihood_mach)
    
    # Fluctuate E[Y|A, W] using a LogisticRegression
    # on the covariate and the offset 
    offset = compute_offset(target_expectation_mach, X)
    covariate = compute_covariate(treatment_likelihood_mach, W, tint)
    fluctuator = glm(reshape(covariate, :, 1), y, Bernoulli(); offset=offset)
    
    # Compute the final estimate tmleATE = 1/n ∑ Fluctuator(t=1, W=w) - Fluctuator(t=0, W=w)
    fluct_treatment_true = compute_fluctuation(fluctuator, 
                                               target_expectation_mach, 
                                               treatment_likelihood_mach, 
                                               W, 
                                               ones(Int, n))
    fluct_treatment_false = compute_fluctuation(fluctuator, 
                                                target_expectation_mach, 
                                                treatment_likelihood_mach, 
                                                W, 
                                                zeros(Int, n))

    tmle.estimate = mean(fluct_treatment_true .- fluct_treatment_false)
    
    # Standard error from the influence curve
    observed_fluct = predict_glm(fluctuator, reshape(covariate, n, 1); offset=offset)
    inf_curve = covariate .* (float(y) .- observed_fluct) .+ fluct_treatment_true .- fluct_treatment_false .- tmle.estimate
    
    tmle.stderror = sqrt(var(inf_curve)/n)

    return tmle

end
