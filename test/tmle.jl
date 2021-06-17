using GenesInteraction
using Random
using Test
using Distributions
using MLJ
using StableRNGs

rng = StableRNG(123)

LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels
DecisionTreeClassifier = @load DecisionTreeClassifier pkg=DecisionTree


function categorical_problem(rng;n=100)
    p_w() = 0.3
    pa_given_w(w) = 1 ./ (1 .+ exp.(-0.5w .+ 1))
    py_given_aw(a, w) = 1 ./ (1 .+ exp.(2w .- 3a .+ 1))
    # Sample from dataset
    Unif = Uniform(0, 1)
    w = rand(rng, Unif, n) .< p_w()
    t = rand(rng, Unif, n) .< pa_given_w(w)
    y = rand(rng, Unif, n) .< py_given_aw(t, w)
    # Convert to dataframe to respect the Tables.jl
    # and convert types
    W = (W=convert(Array{Float64}, w),)
    t = categorical(t)
    y = categorical(y)
    # Compute the theoretical ATE
    ATE₁ = py_given_aw(1, 1)*p_w() + (1-p_w())*py_given_aw(1, 0)
    ATE₀ = py_given_aw(0, 1)*p_w() + (1-p_w())*py_given_aw(0, 0)
    ATE = ATE₁ - ATE₀
    
    return t, W, y, ATE
end


@testset "Test target and treatment should be CategoricalVector{Bool}" begin
    ate_estimator = ATEEstimator(
            LogisticClassifier(),
            LogisticClassifier()
            )
    W = (col1=[1, 1, 0], col2=[0, 0, 1])

    t = categorical([false, true, true])
    y =  categorical(["a", "b", "c"])
    @test_throws MethodError fit!(ate_estimator, t, W, y)

    t = categorical([1, 2, 2])
    y =  categorical([false, true, true])
    @test_throws MethodError fit!(ate_estimator, t, W, y)

end


@testset "Testing the various intermediate functions" begin
    # Let's process the different functions in order of call by fit! 
    # and check intermediate and final outputs
    t = categorical([false, true, true, false, true, true, false])
    y = categorical([true, true, false, false, false, true, true])
    W = (W1=categorical([true, true, false, true, true, true, true]),)
    
    # Testing initial encoding
    hot_mach = machine(OneHotEncoder(), (t=t,))
    fit!(hot_mach)

    X = GenesInteraction.combinetotable(hot_mach, t, W)
    @test X == (t__false=[1, 0, 0, 1, 0, 0, 1],
                t__true=[0, 1, 1, 0, 1, 1, 0],
                W1=W.W1)

    
    # The following dummy model will have the predictive distribution
    # [false, true] => [0.42857142857142855 0.5714285714285714]
    # for both the treatment and the target
    dummy_model = @pipeline OneHotEncoder() ConstantClassifier()
    # Testing compute_offset function
    target_expectation_mach = machine(dummy_model, X, y)
    fit!(target_expectation_mach, verbosity=0)
    offset = GenesInteraction.compute_offset(target_expectation_mach, X)
    @test offset isa Vector{Float64}
    @test all(isapprox.(offset, 0.28768, atol=1e-5))

    # Testing compute_covariate function
    # The likelihood function is given bu the previous distribution
    treatment_likelihood_mach = machine(dummy_model, W, t)
    fit!(treatment_likelihood_mach, verbosity=0)
    covariate = GenesInteraction.compute_covariate(treatment_likelihood_mach, W, t)
    @test covariate isa Vector{Float64}
    @test covariate ≈ [-2.33, 1.75, 1.75, -2.33, 1.75, 1.75, -2.33] atol=1e-2

    # Testing compute_fluctuation function
    fluctuator = GenesInteraction.glm(reshape(covariate, :, 1), y, Bernoulli(); offset=offset)
    fluct = GenesInteraction.compute_fluctuation(fluctuator, 
                                                 target_expectation_mach,
                                                 treatment_likelihood_mach,
                                                 hot_mach,
                                                 W, 
                                                 t)
    # Not sure how to test that the operation is correct
    @test fluct isa Vector{Float64}
    @test all(isapprox.(fluct, [0.6644,
                                0.4977,
                                0.4977,
                                0.6644,
                                0.4977,
                                0.4977,
                                0.6644],
                    atol = 1e-4))
end


@testset "Test ATE TMLE fit! asymptotic behavior" begin
    rng = StableRNG(1234)
    target_cond_expectation_estimator = MLJ.Stack(
                        lr=LogisticClassifier(), 
                        knn=KNNClassifier(),
                        metalearner=LogisticClassifier(), 
                        resampling=StratifiedCV(;nfolds=5, shuffle=false))

    treatment_cond_likelihood_estimator = MLJ.Stack(
                        lr=LogisticClassifier(), 
                        knn=KNNClassifier(),
                        metalearner=LogisticClassifier(), 
                        resampling=StratifiedCV(;nfolds=5, shuffle=false))

    ate_estimator = ATEEstimator(
                        target_cond_expectation_estimator,
                        treatment_cond_likelihood_estimator
                        )
    
    abserrors = []
    for n in [100, 1000, 10000, 100000]
        t, W, y, ATE = categorical_problem(rng; n=n)
        fit!(ate_estimator, t, W, y)
        push!(abserrors, abs(ATE-ate_estimator.estimate))
    end
    @test abserrors == sort(abserrors, rev=true)
end