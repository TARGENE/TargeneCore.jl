
###############################################################################
####################         JLD2Saver Callback            ####################
###############################################################################

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

###############################################################################
####################             Utilities                 ####################
###############################################################################

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

function frequency_criterion(data::DataFrame, columns; threshold=0.001)
    n = size(data, 1)
    freqs = combine(
        groupby(data, columns; skipmissing=true), 
        x -> size(x, 1)/n
    )
    expected_nb_combinations = prod(size(unique(skipmissing(data[!, column])), 1) for column in columns)
    if size(freqs, 1) == expected_nb_combinations && all(x > threshold for x in freqs[!, end])
        return true
    end
    return false
end

"""
If the phenotypes are continuous then we just return all phenotypes
"""
phenotypes_passing_freq_threshold(
    data::DataFrame, 
    ::Type{Real},
    genotypes_columns,
    phenotypes_columns; 
    threshold=0.001) = phenotypes_columns

"""
If the phenotypes are binary then we check frequencies
"""
function phenotypes_passing_freq_threshold(
    data::DataFrame, 
    ::Type{Bool}, 
    genotypes_columns, 
    phenotypes_columns; 
    threshold=0.001
    )
    phenotypes_to_go = []
    for phenotype in phenotypes_columns
        if frequency_criterion(data, [genotypes_columns..., phenotype]; threshold=threshold)
            push!(phenotypes_to_go, phenotype)
        end
    end
    return phenotypes_to_go
end


function preprocess(genotypes, confounders, phenotypes, target_type, queries; 
    freq_threhsold=0.001)
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

    # Check if the lowest frequent genotype passes the threshold
    if !frequency_criterion(data, genotypes_columns; threshold=freq_threhsold)
        return nothing, nothing, nothing, nothing
    end

    # Check if the lowest frequent genotype/phenotype passes the threshold
    phenotypes_columns = phenotypes_passing_freq_threshold(
        data, 
        target_type, 
        genotypes_columns, 
        phenotypes_columns;
        threshold=freq_threhsold
    )
    if isempty(phenotypes_columns)
        return nothing, nothing, nothing, nothing
    end
    # Retrive sample_ids per phenotype
    sample_ids = sample_ids_per_phenotype(data, phenotypes_columns)

    # Retrieve T and convert to categorical data
    # The use of the query he is so that the order of the columns in both
    # T and queries are matching which is currently a technical requirement of the TMLE package
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

###############################################################################
####################           Main Function               ####################
###############################################################################

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
        preprocess(genotypes, confounders, phenotypes, target_type, queries;freq_threhsold=parsed_args["min-freq"])
    # If the frequency threshold has not been met, write a dummy file
    if genotypes === nothing
        jldopen(parsed_args["out"], "a") do io
            JLD2.Group(io, "TMLEREPORTS")
        end
        return 0
    end
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
    return 0
end
