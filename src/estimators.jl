###############################################################################
# BUILD TMLE FROM .TOML


function stack_from_config(config::Dict, metalearner)
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


function estimators_from_toml(config::Dict, queries, run_fn::typeof(PhenotypeTMLEEpistasis))
    tmles = Dict()
    queryvals = [x[2] for x in queries]
    isinteraction = length(queryvals[1]) > 1
    # Parse estimator for the propensity score
    metalearner = SKLogisticClassifier(fit_intercept=false)
    if isinteraction
        G = FullCategoricalJoint(stack_from_config(config["G"], metalearner))
    else
        G = stack_from_config(config["G"], metalearner)
    end
    
    # Parse estimator for the outcome regression
    if haskey(config, "Qcont")
        metalearner =  SKLinearRegressor(fit_intercept=false)
        Q̅ = stack_from_config(config["Qcont"], metalearner)
        tmles["continuous"] = TMLEstimator(Q̅, G, queryvals...)
    end

    if haskey(config, "Qcat")
        metalearner = SKLogisticClassifier(fit_intercept=false)
        Q̅ = stack_from_config(config["Qcat"], metalearner)
        tmles["binary"] = TMLEstimator(Q̅, G, queryvals...)
    end

    length(tmles) == 0 && throw(ArgumentError("At least one of (Qcat, Qcont) "*
                                               "should be specified in the TMLE"*
                                               " configuration file"))
    
    return tmles

end