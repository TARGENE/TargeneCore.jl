save_estimands(outprefix, estimands, batchsize::Nothing) = 
    serialize(batch_name(outprefix, 1), Configuration(estimands=estimands))

function save_estimands(outprefix, estimands, batchsize)
    for (batch_id, batch) ∈ enumerate(Iterators.partition(estimands, batchsize))
        batchfilename = batch_name(outprefix, batch_id)
        serialize(batchfilename, Configuration(estimands=batch))
    end
end

default_order(::Union{typeof(ATE), typeof(CM)}) = [1]
default_order(::typeof(IATE)) = [2]

treatment_tuples_from_single_list(treatment_variables, orders) =
    reduce(vcat, collect(Combinatorics.combinations(treatment_variables, order)) for order in orders)

function estimands_from_flat_list(estimands_configs, dataset, variants, outcomes, confounders; 
    extra_treatments=[],
    outcome_extra_covariates=[],
    positivity_constraint=0.,
    verbosity=1,
    )
    estimands = []
    treatment_variables = isempty(extra_treatments) ? Symbol.(variants) : Symbol.(vcat(variants, extra_treatments))
    for estimands_config in estimands_configs
        estimand_constructor = eval(Symbol(estimands_config["type"]))
        verbosity > 0 && @info(string("Generating estimands: "))
        orders = haskey(estimands_config, "orders") ? estimands_config["orders"] : default_order(estimand_constructor)
        for treatments ∈ treatment_tuples_from_single_list(treatment_variables, orders)
            append!(estimands, factorialEstimands(
                estimand_constructor, dataset,
                treatments, outcomes; 
                confounders=confounders, 
                outcome_extra_covariates=outcome_extra_covariates,
                positivity_constraint=positivity_constraint, 
                verbosity=verbosity-1
            ))
        end
    end
    return estimands
end

function loco_gwas(parsed_args)
    verbosity = parsed_args["verbosity"]
    outprefix = parsed_args["out-prefix"]
    batchsize = parsed_args["batch-size"]
    positivity_constraint = parsed_args["positivity-constraint"]

    loco_gwas_config = parsed_args["loco-gwas"]
    bed_prefix = loco_gwas_config["bed-prefix"]
    traits = read_csv_file(loco_gwas_config["traits"])
    pcs = read_csv_file(loco_gwas_config["pcs"])
    config = YAML.load_file(loco_gwas_config["config"])

    # Variables
    verbosity > 0 && @info("Parsing configuration file.")
    # variants_config = config["variants"]
    extra_treatments = haskey(config, "extra_treatments") ? Symbol.(config["extra_treatments"]) : []
    outcome_extra_covariates = haskey(config, "outcome_extra_covariates") ? Symbol.(config["outcome_extra_covariates"]) : []
    extra_confounders = haskey(config, "extra_confounders") ? Symbol.(config["extra_confounders"]) : []
    confounders = all_confounders(pcs, extra_confounders)
    nonoutcomes = Set(vcat(:SAMPLE_ID, extra_confounders, outcome_extra_covariates, extra_treatments))
    outcomes = filter(x -> x ∉ nonoutcomes, Symbol.(names(traits)))

    # Genotypes and final dataset
    verbosity > 0 && @info("Building and writing dataset.")
    chromosome = SnpData(SnpArrays.datadir(bed_prefix))
    chr_array = convert(Matrix{Float16}, chromosome.snparray)
    genotypes = DataFrame(chr_array, chromosome.snp_info."snpid")
    insertcols!(genotypes, 1, :SAMPLE_ID => chromosome.person_info."iid")

    dataset = merge(traits, pcs, genotypes)
    Arrow.write(string(outprefix, ".data.arrow"), dataset)

    Estimands
    config_type = config["type"]
    estimands_configs = config["estimands"]
    estimands = if config_type == "flat"
        estimands_from_flat_list(estimands_configs, dataset, variants_config, outcomes, confounders;
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity)
    else
        throw(ArgumentError(string("Unknown extraction type: ", type, ", use any of: (flat, groups)")))
    end

    save_estimands(outprefix, groups_ordering(estimands), batchsize)

    verbosity > 0 && @info("Done.")

    return 0
end
