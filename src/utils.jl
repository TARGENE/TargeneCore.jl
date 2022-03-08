#####################################################################
##### GENERIC FUNCTIONS                                         #####
#####################################################################

is_binary(y) = sort(unique(y)) == [0, 1]

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


hdf5filename(::TMLE.Query{<:NamedTuple{names,}}) where names = 
    string(join(names, "_"),  ".hdf5")

jlsfilename(::TMLE.Query{<:NamedTuple{names,}}) where names = 
    string(join(names, "_"),  ".jls")

#####################################################################
#####                 PHENOTYPE PARSING                          ####
#####################################################################

phenotypes_from_data(phenotype_datafile::String) = 
    keys(only(CSV.File(phenotype_datafile, limit=1, drop=["eid"])))

phenotypesnames(phenotype_datafile::String, phenotype_listfile::Nothing) =
    phenotypes_from_data(phenotype_datafile::String)

function phenotypesnames(phenotype_datafile::String, phenotype_listfile::String)
    from_data = phenotypes_from_data(phenotype_datafile)
    from_list = CSV.File(phenotype_listfile, header=[:PHENOTYPES], types=Symbol).PHENOTYPES
    return filter(∈(from_list), from_data)
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

function writeresults(jld_filename, mach::Machine, phenotypename, sample_ids; mach_filename=nothing)
    # Some phenotypes names contain slashes that conflict with JLD2
    phenotypename = replace(String(phenotypename), "/" => "_&_")
    # Save essential information
    jldopen(jld_filename, "a+") do io
        group = JLD2.Group(io, phenotypename)
        group["sample_ids"] = sample_ids
        group["queryreports"] = TMLE.queryreports(mach)
    end
    # Optionnally save the machine
    if mach_filename !== nothing
        serializable!(mach)
        open(mach_filename, "a+") do io
            serialize(io, phenotypename => mach)
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