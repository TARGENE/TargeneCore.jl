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


function retrieve_alleles(probabilities, index_to_genotype, threshold=0.9)
    n = size(probabilities)[2]
    # The default value is missing
    t = Vector{Union{String, Missing}}(missing, n)
    for i in 1:n
        # If no allele has been annotated with sufficient confidence
        # the sample is declared as missing for this variant
        sample_gen_index = findfirst(x -> x >= threshold, probabilities[:, i])
        sample_gen_index isa Nothing || (t[i] = index_to_genotype[sample_gen_index])
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
            index_to_genotype = [alleles(v)[1]^2, alleles(v)[1]*alleles(v)[2], alleles(v)[2]^2]
            probabilities = probabilities!(b, v)
            genotypes[!, rsid] = retrieve_alleles(probabilities, index_to_genotype, threshold)
        end
    end
    return genotypes
end