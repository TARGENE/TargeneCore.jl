

"""
A callback to save the TMLE estimation results to disk instead of keeping them in memory
"""
mutable struct JLD2Saver <: TMLE.Callback
    file::String
    save_machines::Bool
end

function TMLE.after_tmle(callback::JLD2Saver, report::TMLEReport, target_id::Int, query_id::Int)
    jldopen(callback.file, "a") do io
        group = haskey(io, "TMLEREPORTS") ? io["TMLEREPORTS"] : JLD2.Group(io, "TMLEREPORTS")
        group[string(target_id, "_", query_id)] = report
    end
end

function TMLE.after_fit(callback::JLD2Saver, mach::Machine, id::Symbol)
    if callback.save_machines
        jldopen(callback.file, "a") do io
            group = haskey(io, "MACHINES") ? io["MACHINES"] : JLD2.Group(io, "MACHINES")
            group[string(id)] = mach
        end
    end
end

function TMLE.finalize(callback::JLD2Saver, estimation_report::NamedTuple)
    jldopen(callback.file, "a") do io
        io["low_propensity_scores"] = estimation_report.low_propensity_scores
    end
    return estimation_report
end

function sample_ids_per_phenotype(data, phenotypes_columns)
    sample_ids = Dict()
    for colname in phenotypes_columns
        sample_ids[colname] =
            dropmissing(data[!, ["SAMPLE_ID", colname]]).SAMPLE_ID
        
    end
    return sample_ids
end

convert_target(x, ::Type{Bool}) = categorical(x)
convert_target(x, ::Type{Real}) = x

function preprocess(genotypes, confounders, phenotypes, target_type, queries)
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
    for name in keys(first(queries).control)
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
        preprocess(genotypes, confounders, phenotypes, target_type, queries)

    # Build estimator
    tmle_config = TOML.parsefile(parsed_args["estimator"])
    tmle = estimator_from_toml(tmle_config, queries, target_type, adaptive_cv=parsed_args["adaptive-cv"])

    # Run TMLE 
    v >= 1 && @info "TMLE Estimation."
    TMLE.fit(tmle, genotypes, confounders, phenotypes;
             verbosity=v-1, 
             cache=false,
             callbacks=[JLD2Saver(parsed_args["out"], parsed_args["save-full"])]
    )

    # Save sample_ids
    jldopen(parsed_args["out"], "a") do io
        io["SAMPLE_IDS"] = sample_ids
    end

    v >= 1 && @info "Done."
end
