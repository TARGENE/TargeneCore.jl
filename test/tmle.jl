using GenesInteraction
using Random
using Test
using Distributions
using DataFrames
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
    w = rand(rng, Unif, n, 1) .< p_w()
    t = rand(rng, Unif, n) .< pa_given_w(w)
    y = rand(rng, Unif, n) .< py_given_aw(t, w)
    # Convert to dataframe to respect the Tables.jl
    # and convert types
    W = table(convert(Array{Int}, w))
    t = categorical(t)
    y = categorical(y)
    # Compute the theoretical ATE
    ATE₁ = py_given_aw(1, 1)*p_w() + (1-p_w())*py_given_aw(1, 0)
    ATE₀ = py_given_aw(0, 1)*p_w() + (1-p_w())*py_given_aw(0, 0)
    ATE = ATE₁ - ATE₀
    
    return t, W, y, ATE
end


t, W, y, ATE = categorical_problem(rng;n=200)


@testset "Test reformat" begin
    # y should have exactly 2 levels
    W = [1 0;
         1 0;
         0 1]
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
    @test X == [0 1 0;
                1 1 0;
                1 0 1]
    @test tint == [0, 1, 1]
    @test t == tnew
    @test y == ynew
    @test W == Wnew

end


@testset "Test compute_offset" begin
    X = hcat(convert(Vector{Int}, t), W)
    model = LogisticClassifier()
    mach = machine(model, X, y)
    fit!(mach)
    offset = GenesInteraction.compute_offset(mach, X)
end


@testset "Test ATE TMLE fit!" begin
    target_cond_expectation_estimator = MLJ.Stack(
                        lr=LogisticClassifier(), 
                        knn=KNNClassifier(),
                        metalearner=LogisticClassifier(), 
                        resampling=StratifiedCV(;nfolds=3, shuffle=false))

    treatment_cond_likelihood_estimator = MLJ.Stack(
                        lr=LogisticClassifier(), 
                        knn=KNNClassifier(),
                        metalearner=LogisticClassifier(), 
                        resampling=StratifiedCV(;nfolds=3, shuffle=false))

    ate_estimator = ATEEstimator(
                        target_cond_expectation_estimator,
                        treatment_cond_likelihood_estimator
                        )

end