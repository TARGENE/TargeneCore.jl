using MLJ, MLJBase


###################################################
###   SuperLearner and Constructors             ###
###################################################

abstract type AbstractSuperLearner <: DeterministicComposite end

"""
    @superlearner(model1, ..., modelk,
                     metalearner::Model=model,
                     crossval=CV(),
                     name::Symbol=:Superlearner,
                     risk::Function=MLJMeasure)

Implements the SuperLearner estimator described in :
van der Laan, Mark J.; Polley, Eric C.; and Hubbard, Alan E., "Super Learner" (July 2007).
U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 222.
https://biostats.bepress.com/ucbbiostat/paper222

# Current scope

The estimated quantity is E[Y|X]. 
Due to the fact that there is no standard way of identifying the mean of a 
Distribution (`predict_mean` is not enforced):
    - Only Continuous and Binary targets are suported
    - For continuous targets, models in the library should be `Deterministic`


# Required Keywords Arguments
    - name::Symbol=AnyName, it will determine the struct name
    - crossval::{CV, StratifiedCV}=ResamplingStrategy
    - metalearner::Model=AnyModelMatchingTheInputAndTarget

# Examples
```julia
using MLJLinearModels: LogisticClassifier
using GenesInteraction

sl = @superlearner(LogisticClassifier(lambda=1), 
                    metalearner=LogisticClassifier(),
                    crossval=StratifiedCV(;nfolds=10, shuffle=false),
                    name=:TestSuperLearner)
mach = machine(sl, X, y)
fit!(mach)
```
"""
macro superlearner(exs...)
    crossval_ex = :()    
    metalearner_ex = :()
    name_ex = :()
    risk_ex = :(risk=nothing)

    existing_names = []
    library_exs = []
    for ex in exs
        # Parse metalearner
        if typeof(ex) == Expr && ex.args[1] == :metalearner
            model = __module__.eval(ex.args[2])
            metatype = typeof(model)
            metalearner_ex = :(metalearner::$(metatype)=$(model))
        # Parse structure name
        elseif typeof(ex) == Expr && ex.args[1] == :name
            name_ex = ex.args[2].value
        # Parse cross validation strategy
        elseif typeof(ex) == Expr && ex.args[1] == :crossval
            crossval = __module__.eval(ex.args[2])
            crossval_type = typeof(crossval)
            crossval_ex = :(crossval::$(crossval_type)=$(crossval))
        elseif typeof(ex) == Expr && ex.args[1] == :risk
            riskfunc = __module__.eval(ex.args[2])
            risktype = typeof(riskfunc)
            risk_ex = :(risk::$(risktype)=$(riskfunc))
        # Parse library
        else
            model = __module__.eval(ex)
            modeltype = typeof(model)
            attr_name = generate_name!(modeltype, existing_names)
            push!(library_exs, :($(attr_name)::$(modeltype)=$(model)))
        end
    end
    
    struct_ex = :(mutable struct $(name_ex) <: AbstractSuperLearner
        $(metalearner_ex)
        $(crossval_ex)
        $(risk_ex)
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


function getfolds(y::AbstractNode, m::Model, n::Int)
    if m.crossval isa StratifiedCV
        folds = node(YY->MLJBase.train_test_pairs(m.crossval, 1:n, YY), y)
    else
        folds = source(MLJBase.train_test_pairs(m.crossval, 1:n))
    end
    folds
end


function trainrows(X::AbstractNode, folds::AbstractNode, nfold)
    node((XX, ff) -> selectrows(XX, ff[nfold][1]), X, folds)
end


function testrows(X::AbstractNode, folds::AbstractNode, nfold)
    node((XX, ff) -> selectrows(XX, ff[nfold][2]), X, folds)
end


function computerisk(riskfunc, ŷ::AbstractNode, ytrue::AbstractNode)
    node((PP, TT) -> riskfunc(PP, int(TT) ? TT isa CategoricalArray : TT), ŷ, ytrue)
end


function library(m::AbstractSuperLearner)
    filter(x-> x ∉ (:crossval, :metalearner, :risk), fieldnames(typeof(m)))
end


function expected_value(X::AbstractNode)
    node(XX->hasproperty(XX, :prob_given_ref) ? XX.prob_given_ref[2] : XX, X)
end


###################################################
###   Fit function                              ###
###################################################


"""
    MLJ.fit(m::SuperLearner, verbosity::Int, X, y)
"""
function MLJ.fit(m::AbstractSuperLearner, verbosity::Int, X, y)
    n = nrows(y)
    risks = Dict()
    machine_registry = Dict()

    X = source(X)
    y = source(y)

    Zval = []
    yval = []
    
    folds = getfolds(y, m, n)
    for nfold in 1:m.crossval.nfolds
        Xtrain = trainrows(X, folds, nfold)
        ytrain = trainrows(y, folds, nfold)
        Xtest = testrows(X, folds, nfold)
        ytest = testrows(y, folds, nfold)

        Zfold = []
        for modelname in library(m)
            model = getfield(m, modelname)
            mach = machine(model, Xtrain, ytrain)
            machine_registry[(modelname, nfold)] = mach
            ypred = expected_value(predict(mach, Xtest))
            push!(Zfold, ypred)

            if m.risk !== nothing
                risks[(modelname, nfold)] = @node m.risk(ypred, ytest)
            end
        end

        Zfold = hcat(Zfold...)
        
        push!(Zval, Zfold)
        push!(yval, ytest)
    end

    Zval = MLJ.table(vcat(Zval...))
    yval = vcat(yval...)

    metamach = machine(m.metalearner, Zval, yval)
    machine_registry[:metalearner] = metamach

    Zpred = []
    for modelname in library(m)
        model = getfield(m, modelname)
        mach = machine(model, X, y)
        machine_registry[modelname] = mach
        push!(Zpred, expected_value(predict(mach, X)))
    end

    Zpred = MLJ.table(hcat(Zpred...))
    ŷ = expected_value(predict(metamach, Zpred))

    mach = machine(Deterministic(), X, y; predict=ŷ)

    fitresult, cache, report = return!(mach, m, verbosity)

    report = (report..., machine_registry=machine_registry, risks=risks)

    return fitresult, cache, report

end


