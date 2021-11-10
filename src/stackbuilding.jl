###############################################################################
# BUILD TMLE FROM .TOML


function stack_from_config(config::Dict, y)
    # Define the metalearner
    metalearner =  is_binary(y) ? 
        LogisticClassifier(fit_intercept=false) : LinearRegressor(fit_intercept=false)

    # Define the resampling strategy
    resampling = config["resampling"]
    resampling = eval(Symbol(resampling["type"]))(nfolds=resampling["nfolds"])

    # Define the models library
    models = Dict()
    for (modelname, hyperparams) in config
        if !(modelname in ("resampling", "outcome"))
            modeltype = eval(Symbol(modelname))
            paramnames = Tuple(Symbol(x[1]) for x in hyperparams)
            counter = 1
            for paramvals in Base.Iterators.product(values(hyperparams)...)
                model = modeltype(;NamedTuple{paramnames}(paramvals)...)
                models[Symbol(modelname*"_$counter")] = model
                counter += 1
            end
        end
    end

    # Define the Stack
    Stack(;metalearner=metalearner, resampling=resampling, models...)
end


function tmle_from_toml(config::Dict, y)

    F = is_binary(y) ? binaryfluctuation() : continuousfluctuation()
    Q̅ = stack_from_config(config["Q"], y)
    # For now the Treatment is always categorical only
    G = FullCategoricalJoint(stack_from_config(config["G"], [0, 1]))

    return TMLEstimator(Q̅, G, F)

end