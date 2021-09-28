function parse_queries(queryfile::String)
    config = TOML.parsefile(queryfile)
    queries = Dict()
    for (queryname, querydict) in config
        if lowercase(queryname) ∉ ("threshold", "snps")
            rsids = collect(keys(querydict))
            vals = [split(filter(x->!isspace(x), querydict[rsid]), "->") for rsid in rsids]
            rsids_symbols = Tuple(Symbol(x) for x in rsids)
            queries[queryname] = NamedTuple{rsids_symbols}(vals)
        end
    end
    return queries
end


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


function UKBBGenotypes(queryfile)
    config = TOML.parsefile(queryfile)
    snps = config["SNPS"]
    threshold = config["threshold"]
    # Let's load the variants by the files they are in
    bgen_groups = Dict()
    for (rsid, path) in snps
        haskey(bgen_groups, path) ? push!(bgen_groups[path], rsid) : bgen_groups[path] = [rsid]
    end

    genotypes = DataFrame()
    for (path, rsids) in bgen_groups
        b = GenesInteraction.read_bgen(path)

        # Iterate over variants
        for rsid in rsids
            v = variant_by_rsid(b, rsid)
            variant_genotypes = [alleles(v)[1]^2, alleles(v)[1]*alleles(v)[2], alleles(v)[2]^2]
            probabilities = probabilities!(b, v)
            genotypes[!, rsid] = samples_genotype(probabilities, variant_genotypes, threshold)
        end
    end
    return genotypes
end



function prepareTarget(y, config)
    if config["Q"]["outcome"]["type"] == "categorical"
        y = categorical(y)
    end
    return y
end


function filter_data(T, W, y)
    # Combine all elements together
    fulldata = hcat(T, W)
    fulldata[!, "Y"] = y

    # Filter based on missingness
    filtered_data = dropmissing(fulldata)

    return filtered_data[!, names(T)], filtered_data[!, names(W)], filtered_data[!, "Y"]
end

"""
TODO
"""
function UKBBtmleepistasis(phenotypefile, 
    confoundersfile, 
    queryfile, 
    estimatorfile,
    outfile;
    threshold=0.9,
    verbosity=1)
    
    # Build tmle
    tmle = tmle_from_toml(TOML.parsefile(estimatorfile))

    # Parse queries
    queries = parse_queries(queryfile)

    # Build Genotypes
    T = UKBBGenotypes(queryfile; threshold=threshold)

    # Read Confounders
    W = CSV.File(confoundersfile) |> DataFrame

    # Read Target
    y =  DataFrame(CSV.File(phenotypefile;header=false))[:, 3]

    # Filter data based on missingness
    T, W, y = filter_data(T, W, y)

    # Run TMLE over potential epistatic SNPS
    mach = machine(tmle, T, W, y)
    fit!(mach, verbosity)
    println(mach.fitresult.estimate)
    println(mach.fitresult.stderror)
end