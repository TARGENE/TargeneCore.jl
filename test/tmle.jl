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


@testset "Test reformat" begin
    # y should have exactly 2 levels
    W = (col1=[1, 1, 0], col2=[0, 0, 1])

    t = categorical([false, true, true])
    y =  categorical(["a", "b", "c"])
    @test_throws ArgumentError GenesInteraction.reformat(t, W, y)

    # t should be Boolean
    t = categorical([1, 2, 2])
    y =  categorical(["a", "b", "b"])
    @test_throws MethodError GenesInteraction.reformat(t, W, y)

    # Check reformating
    t = categorical([false, true, true])
    y =  categorical(["a", "b", "b"])
    X, tint, tnew, Wnew, ynew = GenesInteraction.reformat(t, W, y)
    @test X == (treatment_ = [0, 1, 1],
                col1 = [1, 1, 0],
                col2 = [0, 0, 1],)
    @test tint == [0, 1, 1]
    @test t == tnew
    @test y == ynew
    @test W == Wnew

end


@testset "Test intermediate computations" begin
    n = 100
    t, W, y, _ = categorical_problem(rng; n=n)
    X, tint, t, W, y = GenesInteraction.reformat(t, W, y)

    # Testing compute_offset function
    target_expectation_mach = machine(LogisticClassifier(), X, y)
    fit!(target_expectation_mach, verbosity=0)
    offset = GenesInteraction.compute_offset(target_expectation_mach, X)
    ## Not sure how to test that the operation is correct
    @test offset isa Vector{Float64}

    # Testing compute_covariate function
    treatment_likelihood_mach = machine(LogisticClassifier(), W, t)
    fit!(treatment_likelihood_mach, verbosity=0)
    covariate = GenesInteraction.compute_covariate(treatment_likelihood_mach, W, tint)
    ## Not sure how to test that the operation is correct
    @test covariate isa Vector{Float64}

    # Testing compute_fluctuation function
    fluctuator = GenesInteraction.glm(reshape(covariate, n, 1), y, Bernoulli(); offset=offset)
    fluct = GenesInteraction.compute_fluctuation(fluctuator, 
                                                 target_expectation_mach,
                                                 treatment_likelihood_mach,
                                                 W, 
                                                 tint)
    ## Not sure how to test that the operation is correct
    @test covariate isa Vector{Float64}
end


@testset "Test ATE TMLE fit!" begin
    n = 10000
    rng = StableRNG(1234)
    t, W, y, ATE = categorical_problem(rng; n=n)
    target_cond_expectation_estimator = MLJ.Stack(
                        lr=LogisticClassifier(), 
                        knn=KNNClassifier(),
                        metalearner=LogisticClassifier(), 
                        resampling=StratifiedCV(;nfolds=10, shuffle=false))

    treatment_cond_likelihood_estimator = MLJ.Stack(
                        lr=LogisticClassifier(), 
                        knn=KNNClassifier(),
                        metalearner=LogisticClassifier(), 
                        resampling=StratifiedCV(;nfolds=10, shuffle=false))

    ate_estimator = ATEEstimator(
                        target_cond_expectation_estimator,
                        treatment_cond_likelihood_estimator
                        )
    
    fit!(ate_estimator, t, W, y)
    # In the large number we can get arbitrarily close
    @test abs(ATE - ate_estimator.estimate) < 0.006

end