module TestInteractionATE
using Test
using GenesInteraction
using MLJ
using Distributions
using Random
using StableRNGs


LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity=0
LinearRegressor = @load LinearRegressor pkg=MLJLinearModels verbosity = 0

function categorical_problem(rng;n)
    p_w() = 0.4
    pt1_given_w(w) = 1 ./ (1 .+ exp.(-0.5w .+ 1))
    pt2_given_w(w) = 1 ./ (1 .+ exp.(0.5w .+ 1))
    py_given_t1t2w(t1, t2, w) = 1 ./ (1 .+ exp.(2w .- 3t1 .+ 3t2 .+ 1))
    # Sample from dataset
    Unif = Uniform(0, 1)
    w = rand(rng, Unif, n) .< p_w()
    t₁ = rand(rng, Unif, n) .< pt1_given_w(w)
    t₂ = rand(rng, Unif, n) .< pt2_given_w(w)
    y = rand(rng, Unif, n) .< py_given_t1t2w(t₁, t₂, w)
    # Format dataset and convert types
    W = (W=convert(Array{Float64}, w),)
    T = (t₁ = categorical(t₁), t₂ = categorical(t₂))
    y = categorical(y)
    # Compute the theoretical ATE
    ATE₁ = (py_given_t1t2w(1, 1, 1) - py_given_t1t2w(1, 0, 1) - py_given_t1t2w(0, 1, 1) + py_given_t1t2w(0, 0, 1))*p_w()
    ATE₀ = (py_given_t1t2w(1, 1, 0) - py_given_t1t2w(1, 0, 0) - py_given_t1t2w(0, 1, 0) + py_given_t1t2w(0, 0, 0))*(1 - p_w())
    ATE = ATE₁ + ATE₀
    
    return T, W, y, ATE
end

@testset "Test tomultivariate" begin
    T = (t1 = categorical([true, false, false, true]),
         t2 = categorical([true, true, false, false]))
    # The mapping is fixed and harcoded
    @test GenesInteraction.tomultivariate(T) == categorical([1, 3, 4, 2])
end

@testset "Testing the various intermediate functions" begin
    T = (t1 = categorical([true, false, false, true, true, true, false]),
         t2 = categorical([true, true, false, false, false, true, false]))
    y = categorical([true, true, false, false, false, true, true])
    W = (W1=categorical([true, true, false, true, true, true, true]),)

    dummy_model = @pipeline OneHotEncoder() ConstantClassifier()
end


@testset "Interaction ATE Asymptotic Behavior" begin

    interaction_estimator = InteractionATEEstimator(
        LogisticClassifier(),
        LogisticClassifier(),
        Bernoulli()
        )

    abs_mean_errors = []
    abs_var_errors = []
    for n in [100, 1000, 10000, 100000]
        abserrors_at_n = []
        for i in 1:10
            rng = StableRNG(i)
            T, W, y, ATE = categorical_problem(rng; n=n)
            fitresult, _, _ = MLJ.fit(interaction_estimator, 0, T, W, y)
            push!(abserrors_at_n, abs(ATE-fitresult.estimate))
        end
        push!(abs_mean_errors, mean(abserrors_at_n))
        push!(abs_var_errors, var(abserrors_at_n))
    end
    # Check the average and variances decrease with n 
    @test abs_mean_errors == sort(abs_mean_errors, rev=true)
    @test abs_var_errors == sort(abs_var_errors, rev=true)
    # Check the error's close to the target for large samples
    @test all(abs_mean_errors .< [0.15, 0.03, 0.02, 0.004])
    @test all(abs_var_errors .< [0.03, 0.0007, 0.0002, 7.9e-6])
end

end

true