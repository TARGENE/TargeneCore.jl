#####################################################################
#####                     QUERY PARSING                          ####
#####################################################################

function parse_queries(queryfile::String)
    config = TOML.parsefile(queryfile)
    queries = Pair[]
    for (queryname, querydict) in config
        if lowercase(queryname) ∉ ("threshold", "snps")
            rsids = collect(keys(querydict))
            vals = [split(filter(x->!isspace(x), querydict[rsid]), "->") for rsid in rsids]
            rsids_symbols = Tuple(Symbol(x) for x in rsids)
            push!(queries, queryname => NamedTuple{rsids_symbols}(vals))
        end
    end
    return queries
end

function querystring(query)
    querystring_ = ""
    for (key, val) in pairs(query)
        querystring_ *= string(key, ": ", join(val, " -> "),  " & ")
    end
    return querystring_[1:length(querystring_) - 3]
end

#####################################################################
#####                 PHENOTYPE PARSING                          ####
#####################################################################

phenotypesnames(phenotypefile::String) = 
    keys(only(CSV.File(phenotypefile, limit=1, drop=["eid"])))


phenotypes_list(phenotype_listfile::Nothing, done_phenotypes, allnames) = 
    filter(x -> x ∉ done_phenotypes, allnames)

function phenotypes_list(phenotype_listfile::String, done_phenotypes, allnames)
    phenotypes_list = CSV.File(phenotype_listfile, 
                               header=[:PHENOTYPES], 
                               types=Symbol)
    filter(x -> (x ∈ Set(phenotypes_list.PHENOTYPES)) & (x ∉ done_phenotypes), allnames)
end



#####################################################################
#####                 CV ADAPTIVE FOLDS                          ####
#####################################################################


"""
    set_cv_folds!(tmle, target, learner; adaptive_cv=true)

Implements the rule of thum given here: https://www.youtube.com/watch?v=WYnjja8DKPg&t=4s
If adaptive_cv is true, it will update the number of folds in the 
cross validation scheme based on the target variable.
"""
function set_cv_folds!(tmle, target; learner=:Q̅, adaptive_cv=false)
    # If adaptive cross validation shouldn't be used
    adaptive_cv == false && return tmle

    # Compute n-eff
    n = nrows(target)
    neff = 
        if scitype(target[1]) == MLJ.Continuous
            n
        else
            counts = [count(==(cat), target) for cat in levels(target)]
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

    # Set the new number of folds
    stack = getfield(tmle, learner)
    stack.resampling = typeof(stack.resampling)(
        nfolds=nfolds,
        shuffle=stack.resampling.shuffle,
        rng=stack.resampling.rng
        )
    
    return tmle
end