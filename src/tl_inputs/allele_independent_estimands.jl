retrieve_variants_list(config::Dict) = vcat((retrieve_variants_list(x) for x in values(config))...)
retrieve_variants_list(variants::AbstractVector) = variants

save_estimands(outprefix, estimands, batchsize::Nothing) = 
    serialize(batch_name(outprefix, 1), Configuration(estimands=estimands))

function save_estimands(outprefix, estimands, batchsize)
    for (batch_id, batch) ∈ enumerate(Iterators.partition(estimands, batchsize))
        batchfilename = batch_name(outprefix, batch_id)
        serialize(batchfilename, Configuration(estimands=batch))
    end
end

"""
    treatment_tuples_from_groups(treatments_lists, orders)

Generate treatment combinations for all order in orders. The final list is sorted 
to encourage reuse of nuisance function in downstreatm estimation.
"""
function treatment_tuples_from_groups(treatments_lists, orders)
    treatment_combinations = []
    for order in orders
        for treatments_lists_at_order in Combinatorics.combinations(treatments_lists, order)
            for treatment_comb in Iterators.product(treatments_lists_at_order...)
                push!(treatment_combinations, treatment_comb)
            end
        end
    end
    return sort(treatment_combinations)
end

function estimands_from_groups(estimands_configs, dataset, variants_config, outcomes, confounders; 
    extra_treatments=[],
    outcome_extra_covariates=[],
    positivity_constraint=0.,
    verbosity=1
    )
    estimands = []
    for (groupname, variants_dict) ∈ variants_config
        verbosity > 0 && @info(string("Generating estimands for group: ", groupname))
        for estimands_config in estimands_configs
            estimand_constructor = eval(Symbol(estimands_config["type"]))
            orders = haskey(estimands_config, "orders") ? estimands_config["orders"] : default_order(estimand_constructor)
            treatments_lists = [Symbol.(variant_list) for variant_list in values(variants_dict)]
            isempty(extra_treatments) || push!(treatments_lists, extra_treatments)
            for treatments ∈ treatment_tuples_from_groups(treatments_lists, orders)
                append!(estimands, factorialEstimands(
                    estimand_constructor, dataset, treatments, outcomes; 
                    confounders=confounders, 
                    outcome_extra_covariates=outcome_extra_covariates,
                    positivity_constraint=positivity_constraint, 
                    verbosity=verbosity-1
                ))
            end
        end
    end
    return estimands
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
    bed_file = loco_gwas_config["bed-file"]
    bim_file = loco_gwas_config["bim-file"]
    fam_file = loco_gwas_config["fam-file"]

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
    chromosome = SnpData(replace(bed_file, r"\.bed$" => ""), famnm=fam_file, bimnm=bim_file)
    chr_array = convert(Matrix{Float16}, chromosome.snparray)
    genotypes = DataFrame(chr_array, chromosome.snp_info."snpid")
    insertcols!(genotypes, 1, :SAMPLE_ID => chromosome.person_info."iid")

    dataset = merge(traits, pcs, genotypes)
    Arrow.write(string(outprefix, ".data.arrow"), dataset)

    variants_config = chromosome.snp_info."snpid"

    #Estimands
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

function allele_independent_estimands(parsed_args)
    verbosity = parsed_args["verbosity"]
    outprefix = parsed_args["out-prefix"]
    batchsize = parsed_args["batch-size"]
    positivity_constraint = parsed_args["positivity-constraint"]

    allele_independent_config = parsed_args["allele-independent"]
    call_threshold = allele_independent_config["call-threshold"]
    bgen_prefix = allele_independent_config["bgen-prefix"]
    traits = read_csv_file(allele_independent_config["traits"])
    pcs = read_csv_file(allele_independent_config["pcs"])
    config = YAML.load_file(allele_independent_config["config"])

    # Variables
    verbosity > 0 && @info("Parsing configuration file.")
    variants_config = config["variants"]
    extra_treatments = haskey(config, "extra_treatments") ? Symbol.(config["extra_treatments"]) : []
    outcome_extra_covariates = haskey(config, "outcome_extra_covariates") ? Symbol.(config["outcome_extra_covariates"]) : []
    extra_confounders = haskey(config, "extra_confounders") ? Symbol.(config["extra_confounders"]) : []
    confounders = all_confounders(pcs, extra_confounders)
    nonoutcomes = Set(vcat(:SAMPLE_ID, extra_confounders, outcome_extra_covariates, extra_treatments))
    outcomes = filter(x -> x ∉ nonoutcomes, Symbol.(names(traits)))

    # Genotypes and final dataset
    verbosity > 0 && @info("Building and writing dataset.")
    variants_set = Set(retrieve_variants_list(variants_config))
    genotypes = call_genotypes(bgen_prefix, variants_set, call_threshold)
    dataset = merge(traits, pcs, genotypes)
    Arrow.write(string(outprefix, ".data.arrow"), dataset)

    # Estimands
    config_type = config["type"]
    estimands_configs = config["estimands"]
    estimands = if config_type == "flat"
        estimands_from_flat_list(estimands_configs, dataset, variants_config, outcomes, confounders;
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity)
    elseif config_type == "groups"
        estimands_from_groups(estimands_configs, dataset, variants_config, outcomes, confounders; 
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity
        )
    else
        throw(ArgumentError(string("Unknown extraction type: ", type, ", use any of: (flat, groups)")))
    end

    save_estimands(outprefix, groups_ordering(estimands), batchsize)

    verbosity > 0 && @info("Done.")

    return 0
end