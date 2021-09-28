module GenesInteraction

using DataFrames
using CSV
using SnpArrays
using MLJ
using Distributions
using TMLE
using TOML
using BGEN

###############################################################################
# INCLUDES

include("ukbb.jl")
include("models.jl")
include("stackbuilding.jl")

###############################################################################
# RUN EPISTASIS ESTIMATION


function prepareTarget(y, config)
    if config["Q"]["outcome"]["type"] == "categorical"
        y = categorical(y)
    end
    return y
end


function filter_data(T, W, y)
    # Combine all elements together
    fulldata = hcat(T, W)
    fulldata[!, "Y"] = y

    # Filter based on missingness
    filtered_data = dropmissing(fulldata)

    return filtered_data[!, names(T)], filtered_data[!, names(W)], filtered_data[!, "Y"]
end

"""
TODO
"""
function UKBBtmleepistasis(phenotypefile, 
    confoundersfile, 
    queryfile, 
    estimatorfile,
    outfile;
    threshold=0.9,
    verbosity=1)
    
    # Build tmle
    config = TOML.parsefile(estimatorfile)
    tmle = tmle_from_toml(config)

    # Build Genotypes
    T = UKBBGenotypes(queryfile; threshold=threshold)

    # Read Confounders
    W = CSV.File(confoundersfile) |> DataFrame

    # Read Target
    y =  DataFrame(CSV.File(phenotypefile;header=false))[:, 3]

    # Filter data based on missingness
    T, W, y = filter_data(T, W, y)

    # Run TMLE over potential epistatic SNPS
    mach = machine(tmle, T, W, y)
    fit!(mach, verbosity)
    println(mach.fitresult.estimate)
    println(mach.fitresult.stderror)
end

###############################################################################
# EXPORTS

export UKBBtmleepistasis


end
