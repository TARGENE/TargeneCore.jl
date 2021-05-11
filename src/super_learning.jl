using MLJ, MLJBase


###################################################
###   SuperLearner and Constructors             ###
###################################################

abstract type AbstractSuperLearner <: DeterministicComposite end

"""
    superlearner model1 model2 ... modelk metalearner::Model=model nfolds::Int=10 shuffle::Bool=false

Implements the SuperLearner estimator described in :
van der Laan, Mark J.; Polley, Eric C.; and Hubbard, Alan E., "Super Learner" (July 2007).
U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 222.
https://biostats.bepress.com/ucbbiostat/paper222

Current Limitations that need further developments:
   - `name` and `metalearner` are required kw arguments
   - Y is binary and a stratified cross validation is used
   - The estimated quantity is E[Y|X]


# Examples
```julia
using MLJLinearModels: LogisticClassifier
using GenesInteraction

superlearner = @superlearner LogisticClassifier(lambda=10) LogisticClassifier(lambda=1) metalearner=LogistiClassifier() nfolds=3 shuffle=false
mach = machine(superlearner, X, y)
fit!(mach)
```
"""
macro superlearner(exs...)
    nfolds_ex = :(nfolds::Int=10)
    shuffle_ex = :(shuffle::Bool=false)
    metalearner_ex = :()
    name_ex = :()
    existing_names = []
    library_exs = []
    for ex in exs
        # Parse shuffle
        if typeof(ex) == Expr && ex.args[1] == :nfolds
            nfolds_ex.args[2] = ex.args[2]
        # Parse nfolds
        elseif typeof(ex) == Expr && ex.args[1] == :shuffle
            shuffle_ex.args[2] = ex.args[2]
        # Parse metalearner
        elseif typeof(ex) == Expr && ex.args[1] == :metalearner
            model = __module__.eval(ex.args[2])
            metatype = typeof(model)
            metalearner_ex = :(metalearner::$(metatype)=$(model))
        elseif typeof(ex) == Expr && ex.args[1] == :name
            name_ex = ex.args[2].value
        # Parse library
        else
            model = __module__.eval(ex)
            modeltype = typeof(model)
            attr_name = generate_name!(modeltype, existing_names)
            push!(library_exs, :($(attr_name)::$(modeltype)=$(model)))
        end
    end
    
    struct_ex = :(mutable struct $(name_ex) <: AbstractSuperLearner
        $(shuffle_ex)
        $(nfolds_ex)
        $(metalearner_ex)
        $(library_exs...)
        end)

    __module__.eval(MLJBase.Parameters.with_kw(struct_ex, __module__, false))

    esc(quote
        $(name_ex)()
        end)
end


###################################################
###   Helper functions                          ###
###################################################


expected_value(X::AbstractNode) = node(XX->XX.prob_given_ref[2], X)

train_test_pairs(X::AbstractNode, cv, rows) = node(XX->MLJBase.train_test_pairs(cv, rows, XX), X)

train_restrict(X::AbstractNode, f::AbstractNode, i) = node((XX, ff) -> restrict(XX, ff[i], 1), X, f)

test_restrict(X::AbstractNode, f::AbstractNode, i) = node((XX, ff) -> restrict(XX, ff[i], 2), X, f)


###################################################
###   Fit function                              ###
###################################################


"""
    MLJ.fit(m::SuperLearner, verbosity::Int, X, y)
"""
function MLJ.fit(m::AbstractSuperLearner, verbosity::Int, X, y)
    stratified_cv = StratifiedCV(; nfolds=m.nfolds, shuffle=m.shuffle)
    n, = size(y)

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
        for modelname in library(m)
            model = getfield(m, modelname)
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
    for modelname in library(m)
        model = getfield(m, modelname)
        mach = machine(model, X, y)
        push!(Zpred, expected_value(predict(mach, X)))
    end

    Zpred = MLJ.table(hcat(Zpred...))
    ŷ = expected_value(predict(metamach, Zpred))

    mach = machine(Deterministic(), X, y; predict=ŷ)

    return!(mach, m, verbosity)

end


function library(m::AbstractSuperLearner)
    filter(x-> x ∉ (:nfolds, :shuffle, :metalearner), fieldnames(typeof(m)))
end