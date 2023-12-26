retrieve_variants_list(config::Dict) = vcat((retrieve_variants_list(x) for x in values(config))...)
retrieve_variants_list(variants::AbstractVector) = variants

mutable struct BatchManager
    outprefix::String
    max_batch_size::Union{Int, Nothing}
    current_batch_id::Int
    current_estimands::Vector
    current_batch_size::Int
    batch_id::Int
    BatchManager(outprefix, max_batch_size) = new(outprefix, max_batch_size, 1, [], 0, 1)
end

function save_batch!(batch_saver::BatchManager, groupname)
    batchfilename = batch_name(string(batch_saver.outprefix, ".", groupname), batch_saver.batch_id)
    serialize(batchfilename, Configuration(estimands=batch_saver.current_estimands))
    batch_saver.batch_id += 1
    batch_saver.current_estimands = []
    batch_saver.current_batch_size = 0
end

function allele_independent_estimands(parsed_args)
    outprefix = parsed_args["out-prefix"]
    batch_saver = BatchManager(outprefix, parsed_args["batch-size"])
    call_threshold = parsed_args["call-threshold"]
    bgen_prefix = parsed_args["bgen-prefix"]
    positivity_constraint = parsed_args["positivity-constraint"]
    traits = read_data(parsed_args["traits"])
    pcs = read_data(parsed_args["pcs"])
    config = YAML.load_file(parsed_args["allele-independent"]["config"])

    # Variables
    variants = config["variants"]
    extra_treatments = haskey(config, "extra_treatments") ? Symbol.(config["extra_treatments"]) : []
    outcome_extra_covariates = haskey(config, "outcome_extra_covariates") ? Symbol.(config["outcome_extra_covariates"]) : []
    extra_confounders = haskey(config, "extra_confounders") ? Symbol.(config["extra_confounders"]) : []
    confounders = all_confounders(pcs, extra_confounders)
    nonoutcomes = Set(vcat(:SAMPLE_ID, extra_confounders, outcome_extra_covariates, extra_treatments))
    outcomes = filter(x -> x ∉ nonoutcomes, Symbol.(names(traits)))

    # Genotypes and final dataset
    variants_set = Set(retrieve_variants_list(variants))
    genotypes = call_genotypes(bgen_prefix, variants_set, call_threshold)
    dataset = merge(traits, pcs, genotypes)
    Arrow.write(string(outprefix, ".data.arrow"), dataset)

    # Estimands
    for (groupname, variants_dict) ∈ variants
        for prod ∈ Iterators.product(values(variants_dict)...)
            treatments = vcat(collect(Symbol.(prod)), extra_treatments)
            treatments_levels = TMLE.unique_treatment_values(dataset, treatments)
            freq_table = positivity_constraint !== nothing ? TMLE.frequency_table(dataset, keys(treatments_levels)) : nothing
            for outcome in outcomes
                Ψ = generateIATEs(
                    treatments_levels, outcome; 
                    confounders=confounders, 
                    outcome_extra_covariates=outcome_extra_covariates,
                    freq_table=freq_table,
                    positivity_constraint=positivity_constraint
                )
                ncomponents = length(Ψ.args)
                if ncomponents > 0 
                    push!(batch_saver.current_estimands, Ψ)
                    batch_saver.current_batch_size += ncomponents + 1
                end
                if batch_saver.max_batch_size !== nothing && batch_saver.current_batch_size > batch_saver.max_batch_size
                    save_batch!(batch_saver, groupname)
                end
            end
        end
        # Save at the end of a group
        if batch_saver.current_batch_size > 0
            save_batch!(batch_saver, groupname)
        end
    end

    return 0
end