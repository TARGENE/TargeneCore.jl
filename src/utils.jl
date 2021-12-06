#####################################################################
#####                 STACK COEFFS PARSING                       ####
#####################################################################

function Qstack_coefs(mach)
    models = getfield(mach.model.Q̅, :models)
    coeffs = fitted_params(mach).Q̅.metalearner.coef
    return [m => c for (m,c) in zip(models, coeffs)]
end

repr_Qstack_coefs(mach; digits=2) = join([m => round(c, digits=digits) for (m, c) in Qstack_coefs(mach)], ", ")

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
                               type=Symbol)
    filter(x -> (x ∈ Set(phenotypes_list.PHENOTYPES)) & (x ∉ done_phenotypes), allnames)
end