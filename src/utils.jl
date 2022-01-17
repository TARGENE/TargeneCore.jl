#####################################################################
#####                 STACK COEFFS PARSING                       ####
#####################################################################

function Qstack_coefs(mach)
    models = getfield(mach.model.Q̅, :models)
    coeffs = [p[2] for p in fitted_params(mach).Q̅.metalearner.coefs]
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
                               types=Symbol)
    filter(x -> (x ∈ Set(phenotypes_list.PHENOTYPES)) & (x ∉ done_phenotypes), allnames)
end

function init_or_retrieve_results(outfile, run_fn::typeof(PhenotypeTMLEEpistasis))
    if isfile(outfile)
        df = CSV.File(outfile, select=[:PHENOTYPE], types=Symbol) |> DataFrame
    else
        df = DataFrame(
            PHENOTYPE=Symbol[],
            QUERYNAME=String[],
            QUERYSTRING=String[],
            ESTIMATE=Float64[], 
            PVALUE=Float64[],
            LOWER_BOUND=Float64[],
            UPPER_BOUND=Float64[],
            STD_ERROR=Float64[],
            QSTACK_COEFS=String[],
            NROWS=Int[],
            TCOUNTS=String[]
            )
        CSV.write(outfile, df)
    end
    return Set(df.PHENOTYPE)
end

function init_or_retrieve_results(outfile, run_fn::typeof(PhenotypeCrossValidation))
    if isfile(outfile)
        df = CSV.File(outfile, select=[:PHENOTYPE], types=Symbol) |> DataFrame
    else
        df = DataFrame(
            PHENOTYPE=Symbol[],
            VARIANTS=String[],
            Q_RESULTSTRING=String[],
            G_RESULTSTRING=String[],
            )
        CSV.write(outfile, df)
    end
    return Set(df.PHENOTYPE)
end

#####################################################################
#####                TMLE REPORT PARSING                         ####
#####################################################################

get_query_report(reports::NamedTuple, i) = reports
get_query_report(reports::Vector, i) = reports[i]


#####################################################################
#####                TREATMENT ENCODING                          ####
#####################################################################

function encode(T)
    n, p = size(T)
    if p == 1
        return T[:, 1]
    end

    joint_levels_it = Iterators.product((levels(T[!, name]) for name in names(T))...)
    encoding = Dict(Tuple(jl) => i for (i, jl) in enumerate(joint_levels_it))

    t_target = Vector{Int}(undef, n)
    for i in 1:n
        t_target[i] = encoding[values(T[i, :])]
    end

    return categorical(t_target; levels=collect(values(encoding)))

end

#####################################################################
#####                 TREATMENT COUNTS                           ####
#####################################################################

function TreatmentCountsRepr(T)
    countrepr = join(names(T), " & ") * ": "
    for row in eachrow(combine(groupby(T, names(T)), nrow))
        countrepr *= join(row, " ") * " | "
    end
    return countrepr[1:end-3]
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