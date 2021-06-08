using GenesInteraction
using Random
using Test
using Distributions
using DataFrames
using MLJ
using StatsPlots


LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels
DecisionTreeClassifier = @load DecisionTreeClassifier pkg=DecisionTree


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

function simulate(ate_estimator, range)
    estimates = []
    trueATES = []
    for n in range
        println("Estimation with n points = $n")
        A, W, y, trueATE = categorical_problem(;n=n)
        ate_estimate = ate_estimator(A, W, y)
        push!(estimates, ate_estimate)
        push!(trueATES, trueATE)
    end

    results = DataFrame(npoints=range,
                        trueATE=trueATES,
                        estimate=estimates)
    return results
end

#######################################################
######## The True model is part of the library ########
#######################################################

target_cond_expectation_estimator = @superlearner(
                    LogisticClassifier(), 
                    KNNClassifier(),
                    metalearner=LogisticClassifier(), 
                    name=:QSuperLearner,
                    crossval=StratifiedCV(;nfolds=3, shuffle=false))

treatment_cond_likelihood_estimator = @superlearner(
                    LogisticClassifier(), 
                    KNNClassifier(),
                    metalearner=LogisticClassifier(), 
                    name=:GSuperLearner,
                    crossval=StratifiedCV(;nfolds=3, shuffle=false))

ate_estimator = ATEEstimator(
                    target_cond_expectation_estimator,
                    treatment_cond_likelihood_estimator
                    )

range = [50, 100, 300, 500, 800, 1000, 2000, 5000, 10000, 20000,
                    50000, 75000, 100000, 200000]
results = simulate(ate_estimator, range)
tmle_estimates = [x.estimate for x in results.estimate]
initial_estimates = [x.initial_estimate for x in results.estimate]
std_errors = [x.standard_error for x in results.estimate]

plot(results.npoints, tmle_estimates;
        ribbon=1.96*std_errors./sqrt.(results.npoints), 
        fillalpha=.5, 
        legend=:bottomright, 
        title="ATE estimates",
        xlabel="Dataset size", 
        ylabel="ATE",
        xaxis=:log,
        label="TML Estimate")
plot!(results.npoints, initial_estimates; label="Initial Estimate")
plot!(results.npoints, results.trueATE; label="Truth")

#######################################################
###### The True model is NOT part of the library ######
#######################################################


# TODO