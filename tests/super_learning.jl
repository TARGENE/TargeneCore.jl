using GenesInteraction
using Distributions, MLJ, Test, DataFrames, CategoricalArrays, MLJLinearModels

function simulation_dataset()
    n = 1000
    Age = rand(DiscreteUniform(30, 70), n)
    X_1 = 1 ./ (1 .+ exp.(-0.005*Age .- 0.01)) .< rand(Uniform(0,1), n)
    X_2 = 1 ./ (1 .+ exp.(-0.01*Age)) .< rand(Uniform(0,1), n)
    X = DataFrame(
        locus_1 = categorical(X_1),
        locus_2 = categorical(X_2),
        age = Age
    )

    y = 1 ./ (1 .+ exp.(-0.01*Age - 0.94*X_1 + 0.4*X_2)) .< rand(Uniform(0,1), n)
    X, categorical(y)
end

library = [@pipeline(OneHotEncoder(), 
                    LogisticClassifier(), 
                    name="MyPipeline")]

X, y = simulation_dataset()
metalearner = GenesInteraction.superlearn(library, X, y;nfolds=3, refit=false)