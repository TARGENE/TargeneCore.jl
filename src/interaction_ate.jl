
"""


"""
mutable struct InteractionATEEstimator <: TMLEstimator 
    target_cond_expectation_estimator::MLJ.Supervised
    treatment_cond_likelihood_estimator::MLJ.Supervised
    fluctuation_family::Distribution
end

function tomultivariate(T)
    t = []
    mapping = Dict((true, true)=>1, (true, false)=>2, (false, true)=>3, (false, false)=>4)
    features = Tables.schema(T).names
    for row in Tables.rows(T)
        push!(t, mapping[(row[features[1]], row[features[2]])])
    end
    return categorical(t)
end


###############################################################################
## Offset and covariate
###############################################################################

"""
(2*T1-1) * (2*T2-1) / p(T1,T2 | W)
"""
function compute_covariate(t_likelihood_estimate::Machine, W, T, t_target)
    # tpred is an estimate of a probability distribution
    # we need to extract the observed likelihood out of it
    features = Tables.schema(T).names
    t₁ = convert(Vector{Int}, Tables.getcolumn(T, features[1]))
    t₂ = convert(Vector{Int}, Tables.getcolumn(T, features[2]))
    tpred = MLJ.predict(t_likelihood_estimate, W)
    likelihood = pdf.(tpred, t_target)
    # truncate predictions, is this really necessary/suitable?
    likelihood = min.(0.995, max.(0.005, likelihood))
    return (2t₁ .- 1).*(2t₂ .- 1) ./ likelihood
end


###############################################################################
## Fluctuation
###############################################################################

function compute_fluctuation(fitted_fluctuator::GeneralizedLinearModel, 
                            target_expectation_mach::Machine, 
                            treatment_likelihood_mach::Machine, 
                            hot_mach::Machine,
                            W, 
                            T, 
                            t_target)
    X = combinetotable(hot_mach, T, W)
    offset = compute_offset(target_expectation_mach, X)
    cov = compute_covariate(treatment_likelihood_mach, W, T, t_target)
    return  GLM.predict(fitted_fluctuator, reshape(cov, :, 1); offset=offset)
end


###############################################################################
## Fit
###############################################################################


function MLJ.fit(tmle::InteractionATEEstimator, 
                 verbosity::Int, 
                 T,
                 W, 
                 y::Union{CategoricalVector{Bool}, Vector{<:Real}})
    n = nrows(y)

    # Fit Encoding of the treatment variable
    t_target = tomultivariate(T)
    hot_mach = machine(OneHotEncoder(), T)
    fit!(hot_mach, verbosity=verbosity)

    # Input checks and reformating
    X = combinetotable(hot_mach, T, W)
    
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

    # Compute the final estimate 
    # InteractionATE = 1/n ∑ [ Fluctuator(t₁=1, t₂=1, W=w) - Fluctuator(t₁=1, t₂=0, W=w)
    #                        - Fluctuator(t₁=0, t₂=1, W=w) + Fluctuator(t₁=0, t₂=0, W=w)]
    counterfactual_treatments = [(true, true, 1), 
                                (true, false, -1), 
                                (false, true, -1), 
                                (false, false, 1)]
    fluct = zeros(n)
    features = Tables.schema(T).names
    for (t₁, t₂, sign) in counterfactual_treatments
        counterfactualT = NamedTuple{features}([
                            categorical(repeat([t₁], n), levels=[true, false]), 
                            categorical(repeat([t₂], n), levels=[true, false])
        ])
        counterfactual_t_target = tomultivariate(counterfactualT)
        fluct .+= sign*compute_fluctuation(fluctuator, 
                                target_expectation_mach, 
                                treatment_likelihood_mach, 
                                hot_mach,
                                W, 
                                counterfactualT,
                                counterfactual_t_target)
    end

    estimate = mean(fluct)

    # Standard error from the influence curve
    observed_fluct = GLM.predict(fluctuator, reshape(covariate, n, 1); offset=offset)
    inf_curve = covariate .* (float(y) .- observed_fluct) .+ fluct .- estimate

    fitresult = (
    estimate=estimate,
    stderror=sqrt(var(inf_curve)/n),
    mean_inf_curve=mean(inf_curve),
    target_expectation_mach=target_expectation_mach,
    treatment_likelihood_mach=treatment_likelihood_mach,
    fluctuation=fluctuator
    )

    cache = nothing
    report = nothing
    return fitresult, cache, report
end