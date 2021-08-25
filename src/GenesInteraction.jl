module GenesInteraction

using DataFrames
using CSV
using SnpArrays
using MLJ
using Distributions
using TMLE
using TOML

###############################################################################
# REGRESSORS

RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree verbosity = 0
NuSVR = @load NuSVR pkg=LIBSVM verbosity = 0
LinearRegressor = @load LinearRegressor pkg=MLJLinearModels verbosity = 0
KNNRegressor = @load KNNRegressor pkg=NearestNeighborModels verbosity = 0
NeuralNetworkRegressor = @load NeuralNetworkRegressor pkg=MLJFlux verbosity = 0
EvoTreeRegressor = @load EvoTreeRegressor pkg=EvoTrees verbosity = 0

###############################################################################
# CLASSIFIERS

RandomForestClassifier= @load RandomForestClassifier pkg=DecisionTree verbosity = 0
NuSVC = @load NuSVC pkg=LIBSVM verbosity = 0
LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity = 0
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels verbosity = 0
NeuralNetworkClassifier = @load NeuralNetworkClassifier pkg=MLJFlux verbosity = 0
EvoTreeClassifier = @load EvoTreeClassifier pkg=EvoTrees verbosity = 0

###############################################################################
# BUILD TMLE FROM .TOML


function stack_from_config(config)
    # Define the metalearner
    outcome = pop!(config, "outcome")
    metalearner = outcome["type"] == "categorical" ? 
        LogisticClassifier(fit_intercept=false) : LinearRegressor(fit_intercept=false)

    # Define the resampling strategy
    resampling = pop!(config, "resampling")
    resampling = eval(Symbol(resampling["type"]))(nfolds=resampling["nfolds"])

    # Define the models library
    models = Dict()
    for (modelname, hyperparams) in config
        modeltype = eval(Symbol(modelname))
        paramnames = Tuple(Symbol(x[1]) for x in hyperparams)
        counter = 1
        for paramvals in Base.Iterators.product(values(hyperparams)...)
            model = modeltype(;NamedTuple{paramnames}(paramvals)...)
            models[Symbol(modelname*"_$counter")] = model
            counter += 1
        end
    end

    # Define the Stack
    Stack(;metalearner=metalearner, resampling=resampling, models...)
end

function tmle_from_toml(file::String)

    config = TOML.parsefile(file)

    ytype = config["Q"]["outcome"]["type"]
    distr = nothing
    if ytype == "categorical"
        distr = Bernoulli()
    elseif ytype == "continuous"
        distr = Normal()
    else
        error("The type of y should be either continuous or categorical")
    end

    Qstack = stack_from_config(config["Q"])
    Gstack = stack_from_config(config["G"])

    return InteractionATEEstimator(Qstack, Gstack, distr)

end




"""
TODO
"""
function runepistatis(genfile, phenotypefile, confoundersfile, snpfile)
    # Read SNP data
    snpdata = SnpData(genfile)
    snp_idx = Dict((x => i) for (i,x) in enumerate(snpdata.snp_info.snpid))
    # Read potential epistatic SNPs
    query_snps = CSV.File(snpfile) |> DataFrame
    # Read Confounders
    W = query_snps = CSV.File(confoundersfile) |> DataFrame
    # Read Target
    Y =  DataFrame(CSV.File(phenotypefile;header=false))[:, 3]
    # Run TMLE over potential epistatic SNPS
    for snps_row in eachrow(query_snps)
        snpcols = [snp_idx[snp] for snp in snps_row]
        T = convert(Matrix{Float64}, @view(snpdata.snparray[:, snpcols]))
        println(size(T))
    end
end

###############################################################################
# EXPORTS

export runepistatis


end
