using Distributions
using MLJ
using GLM: glm
using GLM: predict as predict_glm

abstract type TMLEstimator end

abstract type TMLEstimate end


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
struct ATEEstimator <: TMLEstimator 
    target_cond_expectation_estimator::MLJ.Supervised
    treatment_cond_likelihood_estimator::MLJ.Supervised
end


struct ATEEstimate <: TMLEstimate
    estimate
    initial_estimate
    standard_error
    pvalue
end


function (tmle::ATEEstimator)(A, W, y)
    # Functions that will compute covariates used by the fluctuator
    compute_offset(X) = log.(X ./ (1 .- X))
    compute_covariate(T, L) = reshape((2*T .- 1) ./ L, length(T), 1)

    # A will be used both as a treatment and a target
    # So need to manage types
    A = convert(Vector{Float64}, A)
    A_target = categorical(A)
    X = hcat(A, W)
    n = nrows(y)
    
    # Find an initial estimate of E[Y|A, W] by superlearning
    # and compute the offset
    target_expectation_mach = machine(tmle.target_cond_expectation_estimator, X, y)
    fit!(target_expectation_mach)
    target_expectation = predict(target_expectation_mach)
    offset_full = compute_offset(target_expectation)

    # Find an estimate of P(A|W)
    # and compute covariate
    treatment_likelihood_mach = machine(tmle.treatment_cond_likelihood_estimator, W, A_target)
    fit!(treatment_likelihood_mach)
    treatment_likelihood = predict(treatment_likelihood_mach)
    cov_full = compute_covariate(A, treatment_likelihood)

    # Fluctuate E[Y|A, W] using a LogisticRegression
    fluctuator = glm(cov_full, y, Bernoulli(); offset=offset_full)
    fluct_full = predict_glm(fluctuator, cov_full;offset=offset_full)

    # Compute the final estimate tmleATE = 1/n ∑ Fluctuator(A=1, W=w) - Fluctuator(A=0, W=w)
    treatment_true = ones(n)
    target_treatment_true = predict(target_expectation_mach, hcat(treatment_true, W))
    offset = compute_offset(target_treatment_true)
    cov = compute_covariate(treatment_true, treatment_likelihood)
    fluct_treatment_true = predict_glm(fluctuator, cov;offset=offset)

    treatment_false = zeros(n)
    target_treatment_false = predict(target_expectation_mach, hcat(treatment_false, W))
    offset = compute_offset(target_treatment_false)
    cov = compute_covariate(treatment_false, 1 .- treatment_likelihood)
    fluct_treatment_false = predict_glm(fluctuator, cov; offset=offset)

    tmleATE = mean(fluct_treatment_true .- fluct_treatment_false)
    initialATE = mean(target_treatment_true .- target_treatment_false)

    # Influence curve and standard error
    inf_curve = cov_full .* (float(y) .- fluct_full) .+ fluct_treatment_true .- fluct_treatment_false .- tmleATE
    std_error = sqrt(var(inf_curve)/n)
    pvalue = 2*(1 - cdf(Normal(0, 1), abs(tmleATE/std_error)))
    
    return ATEEstimate(tmleATE, initialATE, std_error, pvalue)

end
