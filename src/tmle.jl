using MLJBase
using MLJ
using MLJLinearModels: LogisticClassifier


mergeAW(A,W) = hcat(A, W)

first_covariate(Qbar) = log.(Qbar ./ (1 .- Qbar))

second_covariate(A, Gbar) = (2*A .- 1) ./ Gbar

compute_ate(ate_1, ate_0) = mean(ate_1 .- ate_0)



function compute_ate(Qbar_estimator, Gbar_estimator, A, W, y)
    # A will be used both as a treatment and a target
    # So need to manage types
    A = convert(Vector{Float64}, A)
    A_target = categorical(A)
    X = mergeAW(A, W)

    n = nrows(y)
    
    # Find an initial estimate of Qbar = E[Y|A, W] by superlearning
    # and compute first covariate X₁
    Qbar_mach = machine(Q_estimator, X, y)
    fit!(Qbar_mach)
    Qbar = predict(Qbar_mach)
    X₁ = first_covariate(Qbar)

    # Find an estimate of Gbar = P(A|W)
    # and compute second covariate X₂
    Gbar_mach = machine(G_estimator, W, A_target)
    fit!(Gbar_mach)
    Gbar = predict(Gbar_mach)
    X₂ = second_covariate(A, Gbar)

    # Fluctuate Qbar = E[Y|A, W]
    # PROBLEM WITH THE INTERCEPT
    covariates = hcat(X₁, X₂)
    fluct_mach = machine(LogisticClassifier(;penalty=:none, fit_intercept=false), covariates, y)
    fit!(fluct_mach)

    # Compute the final estimate ATEₙ = 1/n ∑ Fluctuator(A=1, W=w) - Fluctuator(A=0, W=w)
    A1 = ones(n)
    Qbar_A1 = predict(Qbar_mach, mergeAW(A1, W))
    cov₁ = first_covariate(Qbar_A1)
    cov₂ = second_covariate(A1, Gbar)
    covariates = hcat(cov₁, cov₂)
    fluct_A1 = predict(fluct_mach, covariates)

    A0 = zeros(n)
    Qbar_A0 = predict(Qbar_mach, mergeAW(A0, W))
    cov₁ = first_covariate(Qbar_A0)
    Gbar_ = 1 .- Gbar
    cov₂ = second_covariate(A0, Gbar_)
    covariates = hcat(cov₁, cov₂)
    fluct_A0 = predict(fluct_mach, covariates)

    ATEₙ = mean(fluct_A1 .- fluct_A0)

    return ATEₙ

end
