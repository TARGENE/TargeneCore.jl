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
        if e == ArgumentError("No component passed the positivity constraint.")
            verbosity > 0 && @info "Skipped due to positivity constraint: treatments=$(keys(treatments)), outcomes=$outcomes"
        else
            verbosity > 0 && @warn "Unexpected error creating estimand: treatments=$(keys(treatments)), error=$e"
            throw(e)
        end
    else
        verbosity > 1 && @info "Created $(length(Ψ)) estimands for treatments=$(keys(treatments))"
        append!(estimands, Ψ)
    end
end

"""
    treatment_from_variant(variant, dataset)

    Generate a key-value pair (dictionary) for treatment structs.
"""
function treatments_from_variant(variant::String, dataset::DataFrame)
    variant_levels = sort(levels(dataset[!, variant], skipmissing=true))
    return Dict{Symbol, Vector{UInt8}}(Symbol(variant)=>variant_levels)
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
                treatments = treatments_from_variant(variant, dataset)
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
            if vars_needed > length(extra_treatments) && verbosity > 0
                @warn("Not enough extra treatments to build estimands of order $(order). Need $vars_needed but only $(length(extra_treatments)) provided.")
            end

            for combo in Combinatorics.combinations(extra_treatments, vars_needed)
                for variant in variants
                    if Symbol(variant) in combo
                        @info "Skipping self-interaction: $variant in $combo"
                        continue
                    end
                    treatments = treatments_from_variant(variant, dataset)
                    # Add each extra treatment from the combination
                    for t in combo
                        merge!(treatments, treatments_from_variant(string(t), dataset))
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
        variants_groups = Iterators.partition(variants, length(variants) ÷ Threads.nthreads())
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

get_only_file_with_suffix(files, suffix) = files[only(findall(x -> endswith(x, suffix), files))]

function read_bed_chromosome(bedprefix)
    bed_files = files_matching_prefix(bedprefix)
    fam_file = get_only_file_with_suffix(bed_files, "fam")
    bim_file = get_only_file_with_suffix(bed_files, "bim")
    bed_file = get_only_file_with_suffix(bed_files, "bed")[1:end-4]
    return SnpData(bed_file, famnm=fam_file, bimnm=bim_file)
end

function get_genotypes_from_beds(bedprefix, outprefix, traits, pcs, outcomes)
    snpdata = read_bed_chromosome(bedprefix)
    genotypes = DataFrame(convert(Matrix{UInt8}, snpdata.snparray), snpdata.snp_info."snpid"; makeunique=true)
    genotype_ids = names(genotypes)
    genotype_map = Union{UInt8, Missing}[0, missing, 1, 2]
    for col in genotype_ids
        genotypes[!, col] = [genotype_map[x+1] for x in genotypes[!, col]]
    end
    insertcols!(genotypes, 1, :SAMPLE_ID => snpdata.person_info."iid")
    
    sample_ids = intersect(Set(traits.SAMPLE_ID), Set(pcs.SAMPLE_ID))
    mask = in.(genotypes.SAMPLE_ID, Ref(sample_ids))
    sample_genotypes = genotypes[mask, :]

    counts = countmap.(eachcol(select(sample_genotypes, Not(:SAMPLE_ID))))
    mapping_df = DataFrame(
      snpid   = genotype_ids,
      allele1 = snpdata.snp_info.allele1,
      allele2 = snpdata.snp_info.allele2,
      vₘ      = get.(counts, missing,  0),
      v₀      = get.(counts, UInt8(0), 0),
      v₁      = get.(counts, UInt8(1), 0),
      v₂      = get.(counts, UInt8(2), 0),
    )
    mapping_df.n = mapping_df.v₀ .+ mapping_df.v₁ .+ mapping_df.v₂

    # if v₀ < v₂, swap all 0↔2 so that 0 always marks the major homozygote
    for (i, snpid) in enumerate(mapping_df.snpid)
        if mapping_df.v₀[i] < mapping_df.v₂[i]
            col = Symbol(snpid)
            # Flip genotype coding 0<->2 in both full and sample subsets so downstream per-level counts match mapping
            genotypes_col = genotypes[!, col]
            genotypes[!, col] = map(x -> x === UInt8(0) ? UInt8(2) : x === UInt8(2) ? UInt8(0) : x, genotypes_col)
            if size(sample_genotypes, 1) > 0
                sample_genotypes[!, col] = map(x -> x === UInt8(0) ? UInt8(2) : x === UInt8(2) ? UInt8(0) : x, sample_genotypes[!, col])
            end

            # Swap alleles and counts respectively so that v₀ >= v₂ reflects major allele homozygote
            allele2 = mapping_df.allele2[i]
            mapping_df.allele2[i] = mapping_df.allele1[i]
            mapping_df.allele1[i] = allele2

            v₂_old = mapping_df.v₀[i]
            mapping_df.v₀[i] = mapping_df.v₂[i]
            mapping_df.v₂[i] = v₂_old
        end
    end

    if !isempty(outcomes)
        outcome_syms = Symbol.(outcomes)
        trait_names = names(traits)
        # Keep only outcomes present in traits (either Symbol or String form)
        present_outcomes = [o for o in outcome_syms if (o in trait_names) || (string(o) in trait_names)]
        if isempty(present_outcomes)
            @warn "None of the requested outcomes found in traits; skipping per-outcome genotype counts."
        else
            traits_subset = select(traits, vcat(:SAMPLE_ID, present_outcomes)...)
            joined = innerjoin(sample_genotypes, traits_subset, on = :SAMPLE_ID)
            joined_names = names(joined)

            for outcome in present_outcomes
                outcome_str = string(outcome)
                outcome_col = (outcome in trait_names) ? outcome : outcome_str
                outcome_joined_col = (outcome in joined_names) ? outcome : outcome_str
                lvls = unique(skipmissing(traits[!, outcome_col]))
                if length(lvls) != 2
                    @info "Outcome $(outcome) has $(length(lvls)) levels (not binary); skipping per-level counts."
                    continue
                end
                # Compute counts restricted to genotype columns for each level
                for lvl in lvls
                    rows = joined[joined[!, outcome_joined_col] .== lvl, :]
                    lvl_counts = countmap.(eachcol(select(rows, genotype_ids)))
                    safe_outcome = replace(outcome_str, r"\W+" => "_")
                    safe_lvl = replace(string(lvl), r"\W+" => "_")
                    mapping_df[!, "v0_$(safe_outcome)_$(safe_lvl)"] = get.(lvl_counts, UInt8(0), 0)
                    mapping_df[!, "v1_$(safe_outcome)_$(safe_lvl)"] = get.(lvl_counts, UInt8(1), 0)
                    mapping_df[!, "v2_$(safe_outcome)_$(safe_lvl)"] = get.(lvl_counts, UInt8(2), 0)
                end
            end
        end
    end

    mapping_df.MAF = (mapping_df.v₁ .+ (2 .* mapping_df.v₂)) ./ (2 .* (mapping_df.v₀ .+ mapping_df.v₁ .+ mapping_df.v₂))
    CSV.write("$(outprefix).mapping.txt", mapping_df)
    return genotypes, Dict()
end

function make_genotypes(genotype_prefix, config, call_threshold, outprefix, traits, pcs, outcomes)
    genotypes, genotypes_levels = if config["type"] == "gwas"
        get_genotypes_from_beds(genotype_prefix, outprefix, traits, pcs, outcomes)
    else
        variants_set = Set(retrieve_variants_list(config["variants"]))
        call_genotypes(genotype_prefix, variants_set, call_threshold)
    end
    return genotypes, genotypes_levels
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
    genotypes, treatments_levels = make_genotypes(genotypes_prefix, config, call_threshold, outprefix, traits, pcs, outcomes)
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