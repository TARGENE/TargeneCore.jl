#####################################################################
##### GENERIC FUNCTIONS                                         #####
#####################################################################

isbinary(y) = sort(unique(y)) == [0, 1]

#####################################################################
#####                     QUERY PARSING                          ####
#####################################################################

function parse_queries(queryfile::String)
    config = TOML.parsefile(queryfile)
    queries = Query[]
    for (queryname, querydict) in config
        if lowercase(queryname) ∉ ("threshold", "snps")
            rsids = keys(querydict)
            names = Tuple(Symbol(x) for x in rsids)
            rsid_vals = [split(filter(x->!isspace(x), querydict[rsid]), "->") for rsid in rsids]
            case = NamedTuple{names}(rsid_val[2] for rsid_val in rsid_vals)
            control = NamedTuple{names}(rsid_val[1] for rsid_val in rsid_vals)
            push!(queries,  Query(case=case, control=control, name=queryname))
        end
    end
    return queries
end


#####################################################################
#####                 PHENOTYPE PARSING                          ####
#####################################################################

phenotypesnames(phenotype_listfile::Nothing) = nothing

function phenotypesnames(phenotype_listfile::String)
    phenotypeslist = open(readlines, phenotype_listfile)
    return vcat(["SAMPLE_ID"], phenotypeslist)
end

function load_phenotypes(phenotypes_datafile, phenotypes_listfile)
    columns = phenotypesnames(phenotypes_listfile)
    return CSV.File(phenotypes_datafile, select=columns) |> DataFrame
end

#####################################################################
#####                 CV ADAPTIVE FOLDS                          ####
#####################################################################

countuniques(v::AbstractVector) = [count(==(u), v) for u in unique(v)]
countuniques(table) = 
    countuniques([values(x) for x in Tables.namedtupleiterator(table)])

"""
    AdaptiveCV(cv::Union{CV, StratifiedCV})

Implements the rule of thum given here: https://www.youtube.com/watch?v=WYnjja8DKPg&t=4s
"""
mutable struct AdaptiveCV <: MLJBase.ResamplingStrategy
    cv::Union{CV, StratifiedCV}
end


function MLJBase.train_test_pairs(cv::AdaptiveCV, rows, y)
    # Compute n-eff
    n = nrows(y)
    neff = 
        if scitype(first(y)) == MLJBase.Continuous
            n
        else
            counts = countuniques(y)
            nrare = minimum(counts)
            min(n, 5*nrare)
        end

    # Compute number of folds
    nfolds = 
        if neff < 30
            neff
        elseif neff < 500
            20
        elseif neff < 5000
            10
        elseif neff < 10_000
            5
        else
            3
        end
    
    # Update the underlying n_folds
    adapted_cv = typeof(cv.cv)(nfolds=nfolds, shuffle=cv.cv.shuffle, rng=cv.cv.rng)
    
    return MLJBase.train_test_pairs(adapted_cv, rows, y)
end


#####################################################################
#####                         SAVING                             ####
#####################################################################

function writeresults(outfilename, mach, sample_ids; save_full=false)
    # Save essential information
    jldopen(outfilename, "a+") do io
        io["SAMPLE_IDS"] = sample_ids
        if save_full
            serializable!(mach)
            io["MACHINE"] = mach
        else
            io["QUERYREPORTS"] = queryreports(mach)
        end
    end
end


#####################################################################
#####     ALL THIS SHOULD DISAPPEAR WHEN MLJBASE IS READY        ####
#####################################################################

function wipe_data(mach::Machine)
    mach.data = ()
    mach.args = ()
    mach.resampled_data = ()
end

function serializable!(mach::Machine)
    mach.cache = ()
    wipe_data(mach)
end

function serializable!(mach::Machine{FullCategoricalJoint})
    wipe_data(mach)
    mach.cache = ()
    for submach in machines(glb(mach.fitresult.model_fitresult))
        serializable!(submach)
    end
end

function serializable!(mach::Machine{<:Composite})
    mach.cache = (network_model_names=mach.cache.network_model_names, )
    wipe_data(mach)
    for submach in machines(glb(mach))
        serializable!(submach)
    end
end