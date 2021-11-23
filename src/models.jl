###############################################################################
# REGRESSORS

SKLinearRegressor = @load LinearRegressor pkg=ScikitLearn verbosity=0
KNNRegressor = @load KNNRegressor pkg=NearestNeighborModels verbosity=0
LinearRegressor = @load LinearRegressor pkg=MLJLinearModels verbosity=0
XGBoostRegressor = @load XGBoostRegressor pkg=XGBoost verbosity=0

###############################################################################
# CLASSIFIERS

SKLogisticClassifier = @load LogisticClassifier pkg=ScikitLearn verbosity=0
LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity=0
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels verbosity=0
XGBoostClassifier = @load XGBoostClassifier pkg=XGBoost verbosity=0

###############################################################################
# INTERACTIONS

mutable struct InteractionTransformer <: Unsupervised 
    column_pattern
end


function MLJ.fit(model::InteractionTransformer, verbosity::Int, X)
    matching_columns = filter(
        x -> occursin(model.column_pattern, string(x)), 
        Tables.columnnames(X)
        )
    n = length(matching_columns)
    interaction_pairs = Vector{Pair{Symbol, Symbol}}(undef, n*(n-1)÷2)
    index = 1
    for i in 1:n-1
        for j in i+1:n
            interaction_pairs[index] = matching_columns[i] => matching_columns[j]
            index += 1
        end
    end
    fitresult = (interaction_pairs=interaction_pairs, ninter=length(interaction_pairs))
    return fitresult, nothing, nothing
end


function MLJ.transform(model::InteractionTransformer, fitresult, X)
    n = nrows(X)
    interactions = Matrix{Float64}(undef, n, fitresult.ninter)
    index = 1
    for (col₁, col₂) in fitresult.interaction_pairs
        interactions[:, index] = Tables.getcolumn(X, col₁) .* Tables.getcolumn(X, col₂)
        index += 1
    end
    names = Tuple(Symbol(string(col₁, "_", col₂)) for (col₁, col₂) in fitresult.interaction_pairs)
    interactions_ndt = NamedTuple{names}([interactions[:, i] for i in 1:fitresult.ninter])
    return merge(Tables.columntable(X), interactions_ndt)
end


@pipeline InteractionTransformer(r"^rs[0-9]+") SKLogisticClassifier() name=RSInteractionLogisticClassifier

@pipeline InteractionTransformer(r"^rs[0-9]+") SKLinearRegressor() name=RSInteractionLinearRegressor
