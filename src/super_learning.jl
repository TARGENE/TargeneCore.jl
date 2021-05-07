using MLJ, MLJBase
using MLJLinearModels: LogisticClassifier
using DataFrames


mutable struct SuperLearner <: ProbabilisticNetwork
    library
    metalearner
    nfolds::Int
    shuffle::Bool
end

expected_value(X::AbstractNode) = node(XX->XX.prob_given_ref[2], X)

train_test_pairs(X::AbstractNode, cv, rows) = node(XX->MLJBase.train_test_pairs(cv, rows, XX), X)

train_restrict(X::AbstractNode, f::AbstractNode, i) = node((XX, ff) -> restrict(XX, ff[i], 1), X, f)

test_restrict(X::AbstractNode, f::AbstractNode, i) = node((XX, ff) -> restrict(XX, ff[i], 2), X, f)


function MLJ.fit(m::SuperLearner, verbosity::Int, X, y)
    stratified_cv = StratifiedCV(; nfolds=m.nfolds, shuffle=m.shuffle)
    n, = size(y)
    registry = []

    X = source(X)
    y = source(y)

    Zval = []
    yval = []
    folds = train_test_pairs(y, stratified_cv, 1:n)
    for nfold in 1:m.nfolds
        Xtrain = train_restrict(X, folds, nfold)
        ytrain = train_restrict(y, folds, nfold)
        Xtest = test_restrict(X, folds, nfold)
        ytest = test_restrict(y, folds, nfold)

        Zfold = []
        for model in m.library
            mach = machine(model, Xtrain, ytrain)
            push!(Zfold, expected_value(predict(mach, Xtest)))
        end

        Zfold = hcat(Zfold...)
        
        push!(Zval, Zfold)
        push!(yval, ytest)
    end

    Zval = MLJ.table(vcat(Zval...))
    yval = vcat(yval...)

    metamach = machine(m.metalearner, Zval, yval)

    Zpred = []
    for model in m.library
        mach = machine(model, X, y)
        push!(registry, mach)
        push!(Zpred, expected_value(predict(mach, X)))
    end

    Zpred = MLJ.table(hcat(Zpred...))
    ŷ = predict(metamach, Zpred)

    mach = machine(Probabilistic(), X, y; predict=ŷ)
    fit!(mach, verbosity=verbosity - 1)

    return mach(), nothing, registry

end
