using MLJ, MLJBase
using MLJLinearModels: LogisticClassifier
using DataFrames


"""
Implementation of the SuperLearning routine
"""
function superlearn(library, X, y; nfolds=10, refit=false)
    # Define the learning network
    Xs = source(X)
    ys = source(y)
    Z = DataFrame([])
    for model in library
        m = machine(model, Xs, ys)
        Z_val_fold[!, name(model)] = predict(m, Xs)
    end

    y_hat = nothing

    n, = size(y)
    stratified_cv = StratifiedCV(; nfolds=nfolds, shuffle=false)
    nfolds_iter = MLJBase.train_test_pairs(stratified_cv, 1:n, y)


    Z_val = DataFrame([])
    # Unfortunately can't seem to instantiate empty categorical array
    y_val = nothing

    # Loop over the different folds
    for (train, val) in nfolds_iter
        X_val = X[val, :]
        Z_val_fold = DataFrame([])

        # Loop over the different models in the library
        for model in library
            m = machine(model, X, y)
            fit!(m, rows=train)
            Z_val_fold[!, name(model)] = predict(m, X_val)
        end

        Z_val = vcat(Z_val, Z_val_fold)

        if y_val === nothing
            y_val = y[val]
        else
            y_val = vcat(y_val, y[val])
        end
        
    end

    metalearner = machine(LogisticClassifier(), Z_val, y_val)
    fit!(metalearner, rows=1:n)

    if refit
        for model in library
            m = machine(model, X, y)
            fit!(m)
        end
    end

    metalearner
end


mutable struct SuperLearner <: DeterministicNetwork
    library
    metalearner
end


knc = @load DecisionTreeClassifier pkg=DecisionTree
lr = LogisticClassifier()

library = [knc(), lr]

function MLJ.fit(m::SuperLearner, verbosity::Int, X, y)
    X = source(X)
    y = source(y)

    Z = nothing
    for model in m.library
        mach = machine(model, X, y)

        if Z === nothing
            Z = predict(mach, X)
        else
            hcat(Z, predict(mach, X))
        end
    end

    mach = machine(m.metalearner, Z, y)
    ŷ = predict(mach, Z)

    mach = machine(Deterministic(), X, y; predict=ŷ)
    fit!(mach, verbosity=verbosity - 1)
    return mach()

end






sl = superlearner(library)