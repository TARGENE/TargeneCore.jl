using GenesInteraction
using Random
using Test
using Distributions
using DataFrames
using MLJ

LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels


function categorical_problem(;n=100)
    p_w() = 0.3
    pa_given_w(w) = 1 ./ (1 .+ exp.(-0.5w .+ 1))
    py_given_aw(a, w) = 1 ./ (1 .+ exp.(2w .- 3a .+ 1))
    # Sample from dataset
    Unif = Uniform(0, 1)
    w = rand(Unif, n) .< p_w()
    a = rand(Unif, n) .< pa_given_w(w)
    y = rand(Unif, n) .< py_given_aw(a, w)
    # Convert to dataframe to respect the Tables.jl
    # and convert types
    W = DataFrame(W=convert(Array{Float64, 1}, w))
    A = convert(Vector{Float64}, a)
    y = categorical(y)
    # Compute the theoretical ATE
    ATE₁ = py_given_aw(1, 1)*p_w() + (1-p_w())*py_given_aw(1, 0)
    ATE₀ = py_given_aw(0, 1)*p_w() + (1-p_w())*py_given_aw(0, 0)
    ATE = ATE₁ - ATE₀
    
    return A, W, y, ATE
end


superlearner₁ = @superlearner(
                    LogisticClassifier(), 
                    KNNClassifier(),
                    metalearner=LogisticClassifier(), 
                    name=:QSuperLearner,
                    crossval=StratifiedCV(;nfolds=3, shuffle=false))

superlearner₂ = @superlearner(
    LogisticClassifier(), 
    KNNClassifier(),
    metalearner=LogisticClassifier(), 
    name=:GSuperLearner,
    crossval=StratifiedCV(;nfolds=3, shuffle=false))

    
m = ATETargetedLearningEstimator(superlearner₁, superlearner₂)

A, W, y, true_ate = categorical_problem(;n=100)

mach = machine(m, A, W, y)

fit!(mach)

predict(mach)