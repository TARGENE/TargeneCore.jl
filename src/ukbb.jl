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
    col_filter = ("FID", "SAMPLE_ID")
    genotypes_columns = filter(x -> x ∉ col_filter, names(genotypes))
    confounders_columns = filter(x -> x ∉ col_filter, names(confounders))
    phenotypes_columns = filter(x -> x ∉ col_filter, names(phenotypes))

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

    # Parse queries
    queries = parse_queries(parsed_args["queries"])

    # Output filename
    outfilename = hdf5filename(first(queries))

    v >= 1 && @info "Loading data."
    # Build Genotypes
    genotypes = UKBBGenotypes(parsed_args["queries"], queries)

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
    mach = machine(tmle, genotypes, confounders, phenotypes)
    fit!(mach; verbosity=v-1)

    # Save the results
    writeresults(outfilename, mach, sample_ids; save_full=parsed_args["save-full"])

    v >= 1 && @info "Done."
end
