using Plots: display
using GenesInteraction
using Random
using DataFrames
using Distributions
using MLJ
using Plots


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
    # Reshaping and type coercion
    W = reshape(convert(Array{Float64}, w), n, 1)
    t = categorical(t)
    y = categorical(y)
    # Compute the theoretical ATE
    ATE₁ = py_given_aw(1, 1)*p_w() + (1-p_w())*py_given_aw(1, 0)
    ATE₀ = py_given_aw(0, 1)*p_w() + (1-p_w())*py_given_aw(0, 0)
    ATE = ATE₁ - ATE₀
    
    return t, W, y, ATE
end


function asymptotics_behaviour(ate_estimator, range, rng)
    results = DataFrame(n = Int[],
                        truth = Float64[],
                        tmle = Float64[],
                        lowerbound = Float64[],
                        upperbound = Float64[],
                        baseestimate = Float64[],
                        stderror = Float64[])

    for n in range
        println("Estimation with n points = $n")
        # Generate the data
        t, W, y, ATE = categorical_problem(rng;n=n)
        # Get TMLE
        fit!(ate_estimator, t, W, y)
        # Get base estimate for comparison
        X, _, _, _, _ = GenesInteraction.reformat(t, W, y)
        mach = machine(ate_estimator.target_cond_expectation_estimator, X, y)
        fit!(mach)

        X1 = hcat(ones(n), W)
        y1 = predict(mach, X1).prob_given_ref[2]

        X0 = hcat(zeros(n), W)
        y0 = predict(mach, X0).prob_given_ref[2]

        baseestimate = mean(y1 .- y0)
        # Push results
        result = (n,
                  ATE, 
                  ate_estimator.estimate, 
                  ate_estimator.estimate-1.96ate_estimator.stderror,
                  ate_estimator.estimate+1.96ate_estimator.stderror,
                  baseestimate,
                  ate_estimator.stderror)
        push!(results, result)

    end
    return results
end


function plot_results(results)
    myplot = plot(results.n, results.tmle;
        ribbon=1.96results.stderror, 
        fillalpha=.5, 
        legend=:best, 
        title="ATE estimates",
        xlabel="Dataset size", 
        ylabel="ATE",
        label="TMLE")
    plot!(results.n, results.truth; label="Truth")
    plot!(results.n, results.baseestimate; label="Base Estimate")
    display(myplot)
end


#######################################################
######## The True model is part of the library ########
#######################################################
rng = MersenneTwister(123)
range = [50, 100, 300, 500, 800, 1000, 2000, 5000, 10000]

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

for i in 1:10
    rng = MersenneTwister(i)
    results = asymptotics_behaviour(ate_estimator, range, rng)
    plot_results(results)
end


#######################################################
###### The True model is NOT part of the library ######
#######################################################


target_cond_expectation_estimator = MLJ.Stack(
                    lr=KNNClassifier(), 
                    metalearner=LogisticClassifier(), 
                    resampling=StratifiedCV(;nfolds=3, shuffle=false))

treatment_cond_likelihood_estimator = MLJ.Stack(
                    lr=KNNClassifier(),
                    metalearner=LogisticClassifier(), 
                    resampling=StratifiedCV(;nfolds=3, shuffle=false))

ate_estimator = ATEEstimator(
                    target_cond_expectation_estimator,
                    treatment_cond_likelihood_estimator
                    )

for i in 1:10
    rng = MersenneTwister(i)
    results = asymptotics_behaviour(ate_estimator, range, rng)
    plot_results(results)
end