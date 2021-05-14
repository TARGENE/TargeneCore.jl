using GenesInteraction
using Distributions, Random
using DataFrames, CategoricalArrays
using MLJ

Random.seed!(123)


function simulation_dataset(;n=10000, interaction_coef=0, type=:classification)
    Age = rand(DiscreteUniform(30, 70), n)
    X_1 = 1 ./ (1 .+ exp.(-0.005 * Age .- 0.01)) .< rand(Uniform(0, 1), n)
    X_2 = 1 ./ (1 .+ exp.(-0.01 * Age)) .< rand(Uniform(0, 1), n)
    X = DataFrame(
        locus_1=float(X_1),
        locus_2=float(X_2),
        age=float(Age)
    )

    linear_comb = 0.023 * Age - 0.45 * X_1 + 0.98 * X_2 + interaction_coef * X_1 .* X_2

    if type == :classification
        y = categorical(1 ./ (1 .+ exp.(linear_comb)) .< rand(Uniform(0, 1), n))
    elseif type == :regression
        y = linear_comb + randn(n) * 0.1
    end
    X, y
    
end


function evaluate_all(models, X, y, nfolds=5)
    results = Dict()
    for model in models
        mach = machine(model, X, y)
        e = evaluate!(mach,
                     resampling=CV(rng=1234, nfolds=nfolds), 
                     measure=RootMeanSquaredError())
        μ = round(e.measurement[1], sigdigits=5)
        ste = round(std(e.per_fold[1]) / sqrt(nfolds), digits=5)
        results[MLJ.name(model)] = (μ, ste, mach)
    end
    results
end


struct InteractionAdder <: Static
    col1::String
    col2::String
end


function MLJ.transform(t::InteractionAdder, _, X)
    X[!, "interaction"] = X[!, t.col1] .* X[!, t.col2]
    return X
end


# Dataset generation
n = 10000
interaction_strength = 0.5
X, y = simulation_dataset(;n=n, interaction_coef=interaction_strength, type=:regression)


# Build the SuperLearner

potential_models = models(matching(X, y))

LinearRegressor = @load LinearRegressor pkg=MLJLinearModels
DecisionTreeRegressor = @load DecisionTreeRegressor pkg = DecisionTree
SVMRegressor = @load SVMRegressor pkg = ScikitLearn
KNNRegressor = @load KNNRegressor pkg = NearestNeighborModels
DeterministicConstantRegressor = @load DeterministicConstantRegressor pkg=MLJModels
XGBoostRegressor = @load XGBoostRegressor pkg=XGBoost

sl = @superlearner(LinearRegressor(),
                    SVMRegressor(),
                    DecisionTreeRegressor(),
                    KNNRegressor(), 
                    DeterministicConstantRegressor(),
                    XGBoostRegressor(),
                    metalearner=LinearRegressor(),
                    crossval=CV(;nfolds=3, shuffle=false),
                    name=:RegSuperLearner)


# Evaluation

to_evaluate = [sl, LinearRegressor(), SVMRegressor(), DecisionTreeRegressor(), 
                KNNRegressor(), DeterministicConstantRegressor(), XGBoostRegressor()]

results = evaluate_all(to_evaluate, X, y, 5)

sl_mach = results["RegSuperLearner"][3]

fp = fitted_params(sl_mach.report.machine_registry[:metalearner])
library(sl)

# Adding the true model to the library

pipe_sl = @pipeline InteractionAdder("locus_1", "locus_2") sl name=PipeSL

to_evaluate = [sl, pipe_sl, LinearRegressor(), SVMRegressor(), DecisionTreeRegressor(), 
                KNNRegressor(), DeterministicConstantRegressor(), XGBoostRegressor()]

results = evaluate_all(to_evaluate, X, y, 5)

sl_mach = results["PipeSL"][3].report.machines[2]

fp = fitted_params(sl_mach.report.machine_registry[:metalearner])
library(pipe_sl.reg_super_learner)

lr_params = fitted_params(sl_mach.report.machine_registry[:linear_regressor])

true_params = (coefs = [:locus_1 => -0.45, :locus_2 => 0.98, :age => 0.023, :interaction => interaction_strength], intercept = 0,)
