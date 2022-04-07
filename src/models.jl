###############################################################################
# INTERACTIONS LINEAR MODEL

mutable struct InteractionTransformer <: Unsupervised 
    column_pattern
end


function MLJBase.fit(model::InteractionTransformer, verbosity::Int, X)
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


function MLJBase.transform(model::InteractionTransformer, fitresult, X)
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

struct InteractionLMClassifier <: MLJBase.ProbabilisticComposite
    interaction_transformer::InteractionTransformer
    linear_model
end
struct InteractionLMRegressor <: MLJBase.DeterministicComposite
    interaction_transformer::InteractionTransformer
    linear_model
end

"""
Unfortunately it seems the definition of a pipeline in the module poses problems.
I thus resort to the definition of this trivial learning-network while a better API for pipelines
if not available (Ongoing work as of now: https://github.com/JuliaAI/MLJBase.jl/pull/664).
"""
InteractionLM = Union{InteractionLMClassifier, InteractionLMRegressor}

InteractionLMClassifier(;column_pattern="^rs[0-9]+", kwargs...) =
    InteractionLMClassifier(InteractionTransformer(Regex(column_pattern)), LogisticClassifier(;kwargs...))

InteractionLMRegressor(;column_pattern="^rs[0-9]+", kwargs...) =
    InteractionLMRegressor(InteractionTransformer(Regex(column_pattern)), RidgeRegressor(;kwargs...))


function MLJBase.fit(model::InteractionLM, verbosity::Int, X, y)
    Xs = source(X)
    ys = source(y)

    inter_mach = machine(model.interaction_transformer, Xs)
    Xt = MLJBase.transform(inter_mach, Xs)

    pred_mach = machine(model.linear_model, Xt, ys)
    ŷ = MLJBase.predict(pred_mach, Xt)

    mach = machine(supertype(supertype(typeof(model)))(), Xs, ys; predict=ŷ)

	return!(mach, model, verbosity)
end

MLJBase.target_scitype(model::InteractionLMRegressor) = AbstractVector{<:MLJBase.Continuous}
MLJBase.target_scitype(model::InteractionLMClassifier) = AbstractVector{<:Finite}

MLJBase.input_scitype(model::InteractionLM) = Table