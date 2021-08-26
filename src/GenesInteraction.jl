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
    metalearner = config["outcome"]["type"] == "categorical" ? 
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


function tmle_from_toml(config::Dict)

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


###############################################################################
# RUN EPISTASIS ESTIMATION

"""
    filtering_idx(snparray, snpquery, snp_idx)
Returns a list of indices for which people in snparray have a combination 
of alleles given by the snpquery.
"""
function filtering_idx(snparray, snpquery, snp_idx)
    interaction_order = size(snpquery)[1] ÷ 3
    allindices = BitSet[]
    columns = Int64[]
    for i in 1:interaction_order
        basecol_id = 3i - 2
        snpcol = snp_idx[snpquery[basecol_id]]
        snpindices = findall(
            x -> (x == snpquery[basecol_id + 1] || x == snpquery[basecol_id + 2]), 
            @view(snparray[:, snpcol])
            )
        push!(allindices, BitSet(snpindices))
        push!(columns, snpcol)
    end
    
    return collect(intersect(allindices...)), columns
end

"""

Converts the treatment variables to categorical Table
"""
function prepareTreatment(T, snpquery::DataFrameRow)
    interaction_order = size(snpquery)[1] ÷ 3
    treatment_alleles = [snpquery[3i-1] for i in 1:interaction_order]
    return MLJ.table(T .== transpose(treatment_alleles))
end


function prepareTarget(y, config)
    if config["Q"]["outcome"]["type"] == "categorical"
        y = categorical(y)
    end
    return y
end


"""
TODO
"""
function tmleepistasis(genotypefile, 
    phenotypefile, 
    confoundersfile, 
    snpfile, 
    estimatorfile,
    outfile;
    verbosity=0)

    config = TOML.parsefile(estimatorfile)
    # Build tmle
    tmle = tmle_from_toml(config)
    # Read SNP data
    snpdata = SnpData(genotypefile)
    snp_idx = Dict((x => i) for (i,x) in enumerate(snpdata.snp_info.snpid))
    # Read SNP Queries
    snpqueries = CSV.File(snpfile) |> DataFrame
    # Read Confounders
    W = CSV.File(confoundersfile) |> DataFrame
    # Read Target
    y =  DataFrame(CSV.File(phenotypefile;header=false))[:, 3]
    # Run TMLE over potential epistatic SNPS
    estimates = []
    stderrors = []
    for snpquery in eachrow(snpqueries)
        # Filter and prepare the data to retrieve only individuals 
        # corresponding to the queried SNPs alleles and match `fit`
        # expected input types
        indices, columns = filtering_idx(snpdata.snparray, snpquery, snp_idx)
        Wsnp = W[indices, :]
        ysnp = prepareTarget(y[indices], config)
        Tsnp = prepareTreatment(snpdata.snparray[indices, columns], snpquery)
        # Estimate epistasis
        # TODO: replace with machine syntax
        fitresult, _, _ = TMLE.fit(tmle, verbosity, Tsnp, Wsnp, ysnp)
        push!(estimates, fitresult.estimate)
        push!(stderrors, fitresult.stderror)
    end
    # Save the results
    snpqueries[!, "Estimate"] = estimates
    snpqueries[!, "Standard Error"] = stderrors
    CSV.write(outfile, snpqueries)
end

###############################################################################
# EXPORTS

export tmleepistasis


end
