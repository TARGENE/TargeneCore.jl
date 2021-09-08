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

include("genotype_extraction.jl")
include("models.jl")
include("stackbuilding.jl")

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
    queryfile, 
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
    snpqueries = CSV.File(queryfile) |> DataFrame
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
