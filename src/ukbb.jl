
function sample_ids_per_phenotype(data, phenotypes_columns)
    sample_ids = Dict()
    for colname in  phenotypes_columns
        sample_ids[colname] =
            dropmissing(data[!, ["SAMPLE_ID", colname]]).SAMPLE_ID
        
    end
    return sample_ids
end

convert_target(x, ::Type{Bool}) = categorical(x)
convert_target(x, ::Type{Real}) = x

function preprocess(genotypes, confounders, phenotypes, target_type)
    # Make sure data SAMPLE_ID types coincide
    genotypes.SAMPLE_ID = string.(genotypes.SAMPLE_ID)
    confounders.SAMPLE_ID = string.(confounders.SAMPLE_ID)
    phenotypes.SAMPLE_ID = string.(phenotypes.SAMPLE_ID)
    
    # columns names 
    genotypes_columns = filter(!=("SAMPLE_ID"), names(genotypes))
    confounders_columns = filter(!=("SAMPLE_ID"), names(confounders))
    phenotypes_columns = filter(!=("SAMPLE_ID"), names(phenotypes))

    # Join all elements together by SAMPLE_ID
    data = innerjoin(
            innerjoin(genotypes, confounders, on=:SAMPLE_ID),
            phenotypes,
            on=:SAMPLE_ID
            )

    # Drop missing values from genotypes 
    dropmissing!(data, genotypes_columns)

    # Retrive sample_ids per phenotype
    sample_ids = sample_ids_per_phenotype(data, phenotypes_columns)

    # Retrieve T and convert to categorical data
    T = DataFrame()
    for name in genotypes_columns
        T[:, name] = categorical(data[!, name])
    end

    # Retrieve W
    W = data[!, confounders_columns]

    # Retrieve Y and convert to categorical if needed
    Y = DataFrame()
    for name in phenotypes_columns
        Y[:, name] = convert_target(data[!, name], target_type)
    end

    return T, W, Y, sample_ids
end


function UKBBVariantRun(parsed_args)
    v = parsed_args["verbosity"]
    target_type = parsed_args["target-type"] == "Real" ? Real : Bool

    queries = parse_queries(parsed_args["queries"])

    v >= 1 && @info "Loading data."
    # Load Genotypes
    genotypes = CSV.read(
        parsed_args["genotypes"], 
        DataFrame, 
        select=[:SAMPLE_ID, keys(first(queries).control)...]
    )

    # Read Confounders
    confounders = CSV.File(parsed_args["confounders"]) |> DataFrame

    # Load phenotypes
    phenotypes = load_phenotypes(parsed_args["phenotypes"], parsed_args["phenotypes-list"])

    # Preprocessing
    v >= 1 && @info "Preprocessing data."
    genotypes, confounders, phenotypes, sample_ids = 
        preprocess(genotypes, confounders, phenotypes, target_type)

    # Build estimator
    tmle_config = TOML.parsefile(parsed_args["estimator"])
    tmle = estimator_from_toml(tmle_config, queries, target_type, adaptive_cv=parsed_args["adaptive-cv"])

    # Run TMLE 
    v >= 1 && @info "TMLE Estimation."
    mach = machine(tmle, genotypes, confounders, phenotypes, cache=false)
    fit!(mach; verbosity=v-1)

    # Save the results
    writeresults(parsed_args["out"], mach, sample_ids; save_full=parsed_args["save-full"])

    v >= 1 && @info "Done."
end
