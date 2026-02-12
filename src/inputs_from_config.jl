pcnames(pcs) = filter(!=("SAMPLE_ID"), names(pcs))

confounders_from_pcs(pcs, extraW::Nothing) = Symbol.(pcnames(pcs))
confounders_from_pcs(pcs, extraW) = Symbol.(unique(vcat(pcnames(pcs), extraW)))

retrieve_variants_list(config::Dict) = vcat((retrieve_variants_list(x) for x in values(config))...)
retrieve_variants_list(variants::AbstractVector) = variants

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

function try_append_new_estimands!(
    estimands, 
    dataset, 
    estimand_constructor, 
    treatments, 
    outcomes, 
    confounders;
    outcome_extra_covariates=[],
    positivity_constraint=0.,
    verbosity=1
    )
    local Ψ
    try
        Ψ = factorialEstimands(
        estimand_constructor, treatments, outcomes; 
        confounders=confounders, 
        dataset=dataset,
        outcome_extra_covariates=outcome_extra_covariates,
        positivity_constraint=positivity_constraint, 
        verbosity=verbosity-1
    )
    catch e
        if !(e == ArgumentError("No component passed the positivity constraint."))
            throw(e)
        end
    else
        append!(estimands, Ψ)
    end
end

function gwas_sort(treatment::String, dataset::DataFrame)
    counts = combine(groupby(dataset, treatment, skipmissing=true), nrow)
    genotypes = counts[!, treatment]
    
    # If the treatment is not biallelic default to the vanilla level sort
    all(length.(genotypes) .== 2) || return get_treatment_levels(treatment, dataset; gwas=false)
    
    # Identify the heterozygote
    het = only(filter(g->g[1] != g[2], counts[!, treatment]))
    homo1 = string(het[1], het[1])
    homo2 = string(het[2], het[2])
    geno_counts = Dict(zip(genotypes, counts[!, :nrow]))
    
    major, minor = get(geno_counts, homo1, 0) >= get(geno_counts, homo2, 0) ? (homo1, homo2) : (homo2, homo1)
    treatment_levels = filter(g->g in genotypes, [major, het, minor])

    return Dict{Symbol, Vector{Any}}(Symbol(treatment)=>treatment_levels)

end
"""
    get_treatment_levels(treatment, dataset)

    Generate a key-value pair (dictionary) for treatment structs.
"""
function get_treatment_levels(treatment::String, dataset::DataFrame; gwas::Bool=false)
    if !gwas
        treatment_levels = sort(levels(dataset[!, treatment], skipmissing=true))
        return Dict{Symbol, Vector{Any}}(Symbol(treatment) => treatment_levels)
    else
        gwas_sort(treatment, dataset)
    end
end

function estimands_from_variants(
    variants,
    dataset, 
    estimand_constructor, 
    outcomes, 
    confounders;
    extra_treatments=[],
    outcome_extra_covariates=[],
    positivity_constraint=0.,
    orders=[1],
    verbosity=0
    )

    if estimand_constructor != AIE && estimand_constructor != ATE
        throw(ArgumentError("Using GWAS with estimand type $(estimand_constructor) is not supported. Use ATE or AIE."))
    end

    estimands = []
    for order in orders
        if order == 1
            for variant in variants
                treatments = get_treatment_levels(variant, dataset; gwas=true)
                try_append_new_estimands!(
                    estimands, 
                    dataset, 
                    estimand_constructor, 
                    treatments, 
                    outcomes, 
                    confounders;
                    outcome_extra_covariates=outcome_extra_covariates,
                    positivity_constraint=positivity_constraint,
                    verbosity=verbosity
                )
            end
        else
            vars_needed = order - 1 # Since the variant is already a treatment variable 
            if vars_needed > length(extra_treatments)
                throw(ArgumentError("Not enough extra treatments to build estimands of order $(order). Need $vars_needed but only $(length(extra_treatments)) provided."))
            end

            for combo in Combinatorics.combinations(extra_treatments, vars_needed)
                for variant in variants
                    if Symbol(variant) in combo
                        # Skip combinations where the variant is also in the extra treatments
                        continue
                    end
                    treatments = get_treatment_levels(variant, dataset; gwas=true)
                    # Add each extra treatment from the combination
                    for t in combo
                        if string(t) in variants
                            merge!(treatments, get_treatment_levels(string(t), dataset; gwas=true))
                        else
                            merge!(treatments, get_treatment_levels(string(t), dataset))
                        end
                    end
                    try_append_new_estimands!(
                        estimands, 
                        dataset, 
                        estimand_constructor, 
                        treatments, 
                        outcomes, 
                        confounders;
                        outcome_extra_covariates=outcome_extra_covariates,
                        positivity_constraint=positivity_constraint,
                        verbosity=verbosity
                    )
                end
            end
        end
    end
    return estimands
end

"""
The configuration file specifies extra covariates globally. But some estimands may contain a covariate as a treatment. For example 
sex can be used fro a two points interaction study and is considered a treatment in that case. But in other estimands it is used as a covariate.
We make sure it is not both and the logical thing to do is to remove the treatment from the covariates set if it is contained in it.
"""
get_non_treatment_covariates(covariates, treatments) = setdiff(covariates, collect(treatments))

function estimands_from_groups(estimands_configs, dataset, variants_config, outcomes, confounders, treatments_levels; 
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
                # Remove the potentially duplicated treatments from covariates set
                covariates = get_non_treatment_covariates(outcome_extra_covariates, treatments)
                # Extract the treatment levels for the estimand that is goinf to be created
                estimand_treatments_levels = Dict(T => treatments_levels[T] for T ∈ treatments)
                # Create the estimand if it satisfies the desired marginal positivity constraint
                try_append_new_estimands!(
                    estimands, 
                    dataset, 
                    estimand_constructor, 
                    estimand_treatments_levels, 
                    outcomes, 
                    confounders;
                    outcome_extra_covariates=covariates,
                    positivity_constraint=positivity_constraint,
                    verbosity=verbosity-1
                )
            end
        end
    end
    return estimands
end

default_order(::Union{typeof(ATE), typeof(CM)}) = [1]
default_order(::typeof(AIE)) = [2]

treatment_tuples_from_single_list(treatment_variables, orders) =
    reduce(vcat, collect(Combinatorics.combinations(treatment_variables, order)) for order in orders)

function estimands_from_flat_list(estimands_configs, dataset, variants, outcomes, confounders, treatments_levels; 
    extra_treatments=[],
    outcome_extra_covariates=[],
    positivity_constraint=0.,
    verbosity=1,
    )
    estimands = []
    treatment_variables = isempty(extra_treatments) ? Symbol.(variants) : Symbol.(vcat(variants, extra_treatments))
    for estimands_config in estimands_configs
        estimand_constructor = eval(Symbol(estimands_config["type"]))
        verbosity > 0 && @info(string("Generating estimands: ", estimand_constructor))
        orders = haskey(estimands_config, "orders") ? estimands_config["orders"] : default_order(estimand_constructor)
        for treatments ∈ treatment_tuples_from_single_list(treatment_variables, orders)
            # Extract the treatment levels for the estimand that is goinf to be created
            estimand_treatments_levels = Dict(T => treatments_levels[T] for T ∈ treatments)
            # Create the estimand if it satisfies the desired marginal positivity constraint
            try_append_new_estimands!(
                estimands, 
                dataset, 
                estimand_constructor, 
                estimand_treatments_levels, 
                outcomes, 
                confounders;
                outcome_extra_covariates=outcome_extra_covariates,
                positivity_constraint=positivity_constraint,
                verbosity=verbosity-1
            )
        end
    end
    return estimands
end

function estimands_from_gwas(estimands_configs, dataset, variants, outcomes, confounders; 
    extra_treatments=extra_treatments,
    outcome_extra_covariates = [],
    positivity_constraint=0.,
    verbosity=0
    )
    estimands = []
    for estimands_config in estimands_configs
        estimand_constructor = eval(Symbol(estimands_config["type"]))
        orders = haskey(estimands_config, "orders") ? estimands_config["orders"] : default_order(estimand_constructor)
        partition_size = max(1, cld(length(variants), Threads.nthreads()))
        variants_groups = Iterators.partition(variants, partition_size)
        estimands_tasks = map(variants_groups) do variants
            Threads.@spawn estimands_from_variants(variants, dataset, estimand_constructor, outcomes, confounders;
                extra_treatments=extra_treatments,
                outcome_extra_covariates=outcome_extra_covariates,
                positivity_constraint=positivity_constraint,
                orders=orders,
                verbosity=verbosity
            )
        end
        estimands_partitions = fetch.(estimands_tasks)
        push!(estimands, vcat(estimands_partitions...))
    end
    return vcat(estimands...)
end


function inputs_from_config(config_file, genotypes_prefix, traits_file, pcs_file;
    outprefix="final",
    batchsize=nothing,
    call_threshold=0.9,
    positivity_constraint=0.01,
    verbosity=0)

    traits = read_csv_file(traits_file)
    pcs = load_flash_pca_results(pcs_file)
    config = YAML.load_file(config_file)

    # Variables
    verbosity > 0 && @info("Parsing configuration file.")
    extra_treatments = haskey(config, "extra_treatments") ? Symbol.(config["extra_treatments"]) : []
    outcome_extra_covariates = haskey(config, "outcome_extra_covariates") ? Symbol.(config["outcome_extra_covariates"]) : []
    extra_confounders = haskey(config, "extra_confounders") ? Symbol.(config["extra_confounders"]) : []
    confounders = confounders_from_pcs(pcs, extra_confounders)
    nonoutcomes = Set(vcat(:SAMPLE_ID, extra_confounders, outcome_extra_covariates, extra_treatments))
    outcomes = filter(x -> x ∉ nonoutcomes, Symbol.(names(traits)))

    # Genotypes and final dataset
    verbosity > 0 && @info("Building and writing dataset.")
    genotypes, treatments_levels = make_genotypes(genotypes_prefix, config, call_threshold)
    dataset = merge(traits, pcs, genotypes)
    finalise_treatments_levels!(treatments_levels, extra_treatments, dataset)
    Arrow.write(string(outprefix, ".data.arrow"), dataset)

    # Estimands
    config_type = config["type"]
    estimands = if config_type == "flat"
        estimands_from_flat_list(config["estimands"], dataset, config["variants"], outcomes, confounders, treatments_levels;
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity)
    elseif config_type == "groups"
        estimands_from_groups(config["estimands"], dataset, config["variants"], outcomes, confounders, treatments_levels;
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity
        )
    elseif config_type == "gwas"
        variants = filter(!=("SAMPLE_ID"), names(genotypes))
        # Set default gwas estimands config (ATE) if not specified
        if !haskey(config, "estimands")
            config["estimands"] = [Dict("type" => "ATE")]
        end
        estimands_from_gwas(config["estimands"], dataset, variants, outcomes, confounders;
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity
        )
    else
        throw(ArgumentError(string("Unknown extraction type: ", config_type, ", use any of: (flat, groups, gwas)")))
    end

    @assert length(estimands) > 0 "No estimands left, probably due to a too high positivity constraint."

    write_estimation_inputs(outprefix, dataset, TMLE.groups_ordering(estimands); batchsize=batchsize)

    verbosity > 0 && @info("Done.")

    return 0
end