using MLJ, MLJBase, DataFrames


"""
Implementation of the SuperLearning routine
"""
function superlearn(library, X, y; nfolds=10)
    stratified_cv = StratifiedCV(; nfolds=nfolds, shuffle=false)
    nfolds_iter = MLJBase.train_test_pairs(stratified_cv, 1:length(y), y)

    # Y_val and y_val define a new training set that 
    # will be used to fit the SuperLearner coefficients with :
    # Y_val[i, model_name] = predict(model, X_val[i, :])
    Y_val = DataFrame([])
    y_val = categorical([])

    # Loop over the different folds
    for (train, val) in nfolds_iter
        X_val = X[val, :]
        Y_val_fold = DataFrame([])

        # Loop over the different models in the library
        for model in library
            m = machine(model, X, y)
            fit!(m, rows=train)
            Y_val_fold[!, name(model)] = predict(m, X_val)
        end

        
        Y_val = vcat(Y_val, Y_val_fold)
        y_val = vcat(y_val, y[val])
    end

    
end