function read_bgen(bgen_file::String)
    kwargs = Dict{Symbol, Any}(:sample_path => nothing, :idx_path => nothing)
    if bgen_file[end-3:end] == "bgen"
        base = bgen_file[1:end-4]

        samplefile = base * "sample"
        isfile(samplefile) ? kwargs[:sample_path] = samplefile : nothing

        bgifile = bgen_file * ".bgi"
        isfile(bgifile) ? kwargs[:idx_path] = bgifile : nothing
    end
    return Bgen(bgen_file; kwargs...)
end


function samples_genotype(probabilities, variant_genotypes, threshold=0.9)
    n = size(probabilities)[2]
    # The default value is missing
    t = Vector{Union{String, Missing}}(missing, n)
    for i in 1:n
        # If no allele has been annotated with sufficient confidence
        # the sample is declared as missing for this variant
        sample_gen_index = findfirst(x -> x >= threshold, probabilities[:, i])
        sample_gen_index isa Nothing || (t[i] = variant_genotypes[sample_gen_index])
    end
    return t
end


"""

A heterozygous genotype can be specified as (ALLELE₁, ALLELE₂) or (ALLELE₂, ALLELE₁).
Here we align this heterozygous specification on the query and default to 
(ALLELE₁, ALLELE₂) provided in the BGEN file if nothing is specified in the query.
"""
function variant_genotypes(variant::Variant, queries)
    all₁, all₂ = alleles(variant)
    rsid_ = Symbol(variant.rsid)
    # Either (ALLELE₂, ALLELE₁) is provided in the query
    # and we return it as the heterozygous genotype. 
    # Or the other specification will do in all other cases.
    queried_alleles = [(q.case[rsid_], q.control[rsid_]) for q in queries]
    if all₂*all₁ in collect(Iterators.flatten(queried_alleles))
        return [all₁*all₁, all₂*all₁, all₂*all₂]
    end
    return [all₁*all₁, all₁*all₂, all₂*all₂]
end


function UKBBGenotypes(queryfile, queries)
    config = TOML.parsefile(queryfile)
    snps = config["SNPS"]
    threshold = config["threshold"]
    # Let's load the variants by the files they are in
    bgen_groups = Dict()
    for (rsid, path) in snps
        haskey(bgen_groups, path) ? push!(bgen_groups[path], rsid) : bgen_groups[path] = [rsid]
    end

    genotypes = nothing
    for (path, rsids) in bgen_groups
        b = read_bgen(path)
        chr_genotypes = DataFrame(SAMPLE_ID=b.samples)

        # Iterate over variants in this chromosome
        for rsid in rsids
            v = variant_by_rsid(b, rsid)
            variant_gens = variant_genotypes(v, queries)
            probabilities = probabilities!(b, v)
            chr_genotypes[!, rsid] = samples_genotype(probabilities, variant_gens, threshold)
        end
        # I think concatenating should suffice but I still join as a safety
        genotypes isa Nothing ? genotypes = chr_genotypes :
            genotypes = innerjoin(genotypes, chr_genotypes, on=:SAMPLE_ID)
    end
    # Reorder the SNPs to match the order of the queries
    ordered_snps = keys(first(queries).case)
    return genotypes[!, [:SAMPLE_ID, ordered_snps...]]
end


function preprocess(genotypes, confounders, phenotypes;
                    verbosity=1)
    
    # Make sure data SAMPLE_ID types coincide
    genotypes.SAMPLE_ID = string.(genotypes.SAMPLE_ID)
    confounders.SAMPLE_ID = string.(confounders.SAMPLE_ID)
    phenotypes.eid = string.(phenotypes.eid)
    
    # Join all elements together
    data = innerjoin(
            innerjoin(genotypes, confounders, on=:SAMPLE_ID),
            phenotypes,
            on = :SAMPLE_ID =>:eid
            )

    # Check no data has been lost at this stage
    any(nrows(array) != nrows(data) for array in (genotypes, confounders, phenotypes)) &&
        @warn """The number of matching samples is different in data sources: \n 
            - Genotypes: $(nrows(genotypes)) \n
            - Confounders: $(nrows(confounders)) \n
            - Phenotypes: $(nrows(phenotypes)) \n
            - After Join: $(nrows(data)) \n
        """

    # Filter any line where an element is missing
    filtered_data = dropmissing(data)
    verbosity >= 1 && @info "Samples size after missing removed: $(nrows(filtered_data))"

    # Retrieve T and convert to categorical data
    T = filtered_data[!, filter(!=("SAMPLE_ID"), DataFrames.names(genotypes))]
    for name in DataFrames.names(T)
        T[!, name] = categorical(T[:, name])
    end

    # Retrieve W
    W = filtered_data[!, filter(!=("SAMPLE_ID"), names(confounders))]

    # Retrieve y which is assumed to be the first column after SAMPLE_ID
    # and convert to categorical array if needed
    y = filtered_data[!, filter(!=("eid"), names(phenotypes))][!, 1]
    
    # Only support binary and continuous traits for now
    is_binary(y) ? y = categorical(convert(Vector{Bool}, y)) : nothing

    return T, W, y
end


function PhenotypeTMLEEpistasis(tmles::Dict, T, W, y; verbosity=1)
    # Build TMLE machine, can I update Y withoug altering the machine state
    # so that only the Q mach is fit again?
    tmle = is_binary(y) ? tmles["binary"] : tmles["continuous"]
    mach = machine(tmle, T, W, y)
    # Run TMLE 
    fit!(mach; verbosity=verbosity-1)

    return mach
end


function UKBBVariantRun(parsed_args)
    v = parsed_args["verbosity"]

    # Parse queries
    queries = parse_queries(parsed_args["queries"])

    # Build estimators
    tmle_config = TOML.parsefile(parsed_args["estimator"])
    estimators = estimators_from_toml(tmle_config, queries)

    v >= 1 && @info "Loading Genotypes and Confounders."
    # Build Genotypes
    genotypes = UKBBGenotypes(parsed_args["queries"], queries)

    # Read Confounders
    confounders = CSV.File(parsed_args["confounders"]) |> DataFrame

    # Generate phenotypes range
    phenotypenames_ = phenotypesnames(parsed_args["phenotypes"], parsed_args["phenotypes-list"])
    open(parsed_args["output"], "w") do file
        for phenotypename in phenotypenames_
            v >= 1 && @info "Running procedure with phenotype: $phenotypename."
            # Read Target
            phenotype = CSV.File(parsed_args["phenotypes"], select=[:eid, phenotypename]) |> DataFrame
            # Preprocess all data together
            T, W, y = preprocess(genotypes, 
                                confounders, 
                                phenotype;
                                verbosity=v)
            # Update Cross validation settings
            set_cv_folds!(estimators, y, learner=:Q̅, adaptive_cv=parsed_args["adaptive-cv"])
            # Run the estimation procedure
            mach = PhenotypeTMLEEpistasis(estimators, T, W, y; verbosity=v)
            # Update the results
            writeresults(file, mach, phenotypename, full=parsed_args["savefull"])
        end
    end

    v >= 1 && @info "Done."
end
