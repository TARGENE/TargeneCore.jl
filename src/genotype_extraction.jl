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


function retrieve_alleles(probabilities, all_1_index::Int, all_2_index::Int; threshold=0.9)
    n = size(probabilities)[2]
    # The default value is missing
    t = Vector{Union{Bool, Nothing}}(undef, n)
    for i in 1:n
        # If no allele has been annotated with sufficient confidence
        # we discard the sample
        sample_allele_index = findfirst(x -> x > threshold, probabilities[:, i])
        if sample_allele_index == all_1_index
            t[i] = true
        elseif sample_allele_index == all_2_index
            t[i] = false
        end
    end
    return t
end


function variant_alleles_indices(v::Variant, alleles_str::String)
    phased(v) == true && throw("The current assumption is that the data is not phased.")

    allele_counts_to_idx = Dict((2, 0) => 1, (1, 1) => 2, (0, 2) => 3)
    v_alleles = alleles(v)
    allellestr_to_idx = Dict(v_alleles[i] => i for i in 1:size(v_alleles)[1])
    genotype = (allellestr_to_idx[x] for x in split(alleles_str, '/'))
    allele_counts = (sum(x==1 for x in genotype), sum(x==2 for x in genotype))
    return allele_counts_to_idx[allele_counts]
end


function UKBBGenotypes(queryfile; threshold=0.9)
    query_df = CSV.File(queryfile) |> DataFrame
    # Let's load the variants by the files they are in
    # Not sure how useful that is...
    bgen_groups = groupby(query_df, :BGENFILE)
    invalid_ids = Set{Int}()
    genotypes = DataFrame()
    for key in keys(bgen_groups)
        b = GenesInteraction.read_bgen(key.BGENFILE)
        group = bgen_groups[key]
        # Iterate over variants
        for row in eachrow(group)
            v = variant_by_rsid(b, row.RSID)
            probabilities = probabilities!(b, v)
            all_1_index = variant_alleles_indices(v, row.ALLELE_1)
            all_2_index = variant_alleles_indices(v, row.ALLELE_2)
            genotypes[!, row.RSID] = retrieve_alleles(probabilities, all_1_index, all_2_index; threshold=threshold)
        end
    end
    return genotypes
end