const CHR_REG = r"chr[1-9]+"

check_chr_format(x) = @assert all(occursin.(CHR_REG, x))

function allele_specific_binding_snps(asb_prefix)
    bQTLs = DataFrame()
    asbdir_, prefix_ = splitdir(asb_prefix)
    asbdir = asbdir_ == "" ? "." : asbdir_
    for file in readdir(asbdir)
        if occursin(prefix_, file)
            new_bQTLs = CSV.read(joinpath(asbdir_, file), DataFrame, select=[:ID, :CHROM, :isASB])
            bQTLs = vcat(
                bQTLs, 
                DataFrames.select(filter(x -> x.isASB === true, new_bQTLs), [:ID, :CHROM])
            )
        end
    end
    check_chr_format(bQTLs.CHROM)
    bQTLs.QTL_TYPE .= :bQTL
    return unique(bQTLs, :ID)
end

function trans_actors(path)
    eQTLs = CSV.read(path, DataFrame, select=[:ID, :CHROM])
    check_chr_format(eQTLs.CHROM)
    eQTLs.QTL_TYPE .= :eQTL
    return eQTLs
end

function add_chrfiles(snps, chr_prefix)
    chr_dir_, prefix_ = splitdir(chr_prefix)
    chr_dir = chr_dir_ == "" ? "." : chr_dir_
    required_chrs = unique(snps.CHROM)
    chr_files = DataFrame(CHROM=[], BGEN_FILE=[])
    for filename in readdir(chr_dir)
        if occursin(prefix_, filename) && endswith(filename, "bgen")
            chromosome = match(CHR_REG, filename).match
            if chromosome ∈ required_chrs
                filepath = joinpath(chr_dir_, filename)
                sample_filepath = string(filepath[1:end-4], "sample")
                idx_filepath = string(filepath, ".bgi")
                bgenfile = Bgen(filepath, sample_path=sample_filepath, idx_path=idx_filepath)
                push!(chr_files, [chromosome, bgenfile])
            end
        end
    end
    nb_required_chrs = size(required_chrs, 1)
    nb_actual_chrs = size(chr_files, 1)
    @assert nb_required_chrs == nb_actual_chrs "Required "*
        "$nb_required_chrs chromosome files, "* 
        "got: $nb_actual_chrs. Chromosome files: \n $(chr_files.CHROM)"
    return innerjoin(snps, chr_files, on=:CHROM)
end

exclude!(snps, exclusion_file::Nothing) = nothing
exclude!(snps, exclusion_file) = filter!(:ID => x -> x ∉ readlines(exclusion_file), snps)

function impute_genotypes(probabilities, variant_genotypes, threshold)
    n = size(probabilities)[2]
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
Only bi-allelic SNPs are considered. We enforce that the allele₁, allele₍
genotype is encoded with the major variant coming first.
"""
function genotypes(variant::Variant, minor, major)
    alleles_ = alleles(variant)
    @assert size(alleles_, 1) == 2 "$(variant.rsid) is not bi-allelic"
    all₁, all₂ = alleles_
    @assert (all₁, all₂) == (minor, major) || (all₁, all₂) == (major, minor)
    return [all₁*all₁, major*minor, all₂*all₂]
end


function impute_genotypes(snps, threshold)
    snps_genotypes = nothing
    rsid_to_minor_major = Dict()
    for group in groupby(snps, :CHROM)
        bgenfile = first(group.BGEN_FILE)
        chr_genotypes= DataFrame(SAMPLE_ID=bgenfile.samples)
        # Iterate over variants in this chromosome
        for rsid in group.ID
            v = variant_by_rsid(bgenfile, rsid)
            minor_allele_dosage!(bgenfile, v)
            major = major_allele(v)
            minor = minor_allele(v)
            rsid_to_minor_major[rsid] = (minor=minor, major=major)
            variant_genotypes = genotypes(v, minor, major)
            probabilities = probabilities!(bgenfile, v)
            chr_genotypes[!, rsid] = impute_genotypes(probabilities, variant_genotypes, threshold)
        end
        # I think concatenating should suffice but I still join as a safety
        snps_genotypes = snps_genotypes isa Nothing ? chr_genotypes :
            innerjoin(snps_genotypes, chr_genotypes, on=:SAMPLE_ID)
    end
    return snps_genotypes, rsid_to_minor_major
end

function build_queries(eQTLs, bQTLs, genotypes, rsid_to_major_minor, threshold)
    queries = []
    n = size(genotypes, 1)
    for eqtl in eachrow(eQTLs)
        # Control is MAJOR*MAJOR
        # Treatment is MAJOR*MINOR
        eqtl_minor, eqtl_major = rsid_to_major_minor[eqtl.ID]
        eqtl_control = eqtl_major*eqtl_major
        eqtl_case = eqtl_major*eqtl_minor
        for bqtl in eachrow(bQTLs)
            bqtl_minor, bqtl_major = rsid_to_major_minor[bqtl.ID]
            bqtl_control = bqtl_major*bqtl_major
            bqtl_case = bqtl_major*bqtl_minor

            counts = combine(
                groupby(genotypes, [eqtl.ID, bqtl.ID]; skipmissing=true), 
                x -> size(x, 1)
            )
            of_interest = filter(x-> x[1] ∈ (eqtl_control, eqtl_case) && x[2] ∈ (bqtl_control, bqtl_case), counts)
            if size(of_interest, 1) == 4 && minimum(of_interest[!, 3])/n >= threshold
                push!(queries, Dict(
                        "SNPS" => Dict(eqtl.ID => eqtl.CHROM, bqtl.ID => bqtl.CHROM),
                        "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR" => Dict(
                            eqtl.ID => eqtl_control*" -> "*eqtl_case,
                            bqtl.ID => bqtl_control*" -> "*bqtl_case
                        )
                    )
                )
            end
        end
    end
    return queries
end

function write_outputs(genotypes, queries, outdir)
    CSV.write(joinpath(outdir, "genotypes.csv"), genotypes)
    for query in queries
        query_filename = "query_"
        for (snp, _) in query["SNPS"]
            query_filename = string(query_filename, snp, "__")
        end
        query_filename = string(query_filename[1:end-2], ".toml")
        open(joinpath(outdir, query_filename), "w") do io
            TOML.print(io, query)
        end
    end
end

function read_queries(query_prefix)
    queries = Dict[]
    querydir_, prefix = splitdir(query_prefix)
    querydir = querydir_ == "" ? "." : querydir_
    for filename in readdir(querydir)
        if occursin(prefix, filename)
            push!(queries, TOML.parse(open(joinpath(querydir_, filename))))
        end
    end
    return queries
end


function snps_from_queries(queries)
    snps = DataFrame(ID=String[], CHROM=String[])
    for query in queries
        for (snp, chr) in query["SNPS"]
            push!(snps, (snp, chr))
        end
    end
    return unique(snps, :ID)
end


function maybe_filter()

end

function asb_generated_mode(parsed_args)
    # Load SNPS
    eQTLs = trans_actors(parsed_args["trans-actors"])
    bQTLs = allele_specific_binding_snps(parsed_args["asb-prefix"])
    # Build genotypes
    snps = vcat(eQTLs, bQTLs)
    exclude!(snps, parsed_args["exclude"])
    snps = add_chrfiles(snps, parsed_args["chr-prefix"])
    genotypes, rsid_to_minor_major = impute_genotypes(snps, parsed_args["call-threshold"])
    # Build queries
    eQTLs = filter(:QTL_TYPE => ==(:eQTL), snps)
    bQTLs = filter(:QTL_TYPE => ==(:bQTL), snps)
    queries = build_queries(eQTLs, bQTLs, genotypes, rsid_to_minor_major, parsed_args["minor-genotype-freq"])
    # write genotypes and queries
    write_outputs(genotypes, queries, parsed_args["outdir"])
end


function given_queries_mode(parsed_args)
    # Read provided queries
    queries = read_queries(parsed_args["query-prefix"])
    # Get snps from queries and add chromosome file information
    snps = snps_from_queries(queries)
    exclude!(snps, parsed_args["exclude"])
    snps = add_chrfiles(snps, parsed_args["chr-prefix"])
    # Get required genotypes
    genotypes, _ = impute_genotypes(snps, parsed_args["call-threshold"])
    # write genotypes and rewrite queries to outdir
    # maybe in the future check that genotypes in queries 
    #are not too rare as for the asb mode
    write_outputs(genotypes, queries, parsed_args["outdir"])
end


function build_genotypes_and_queries(parsed_args)
    if parsed_args["mode"] == "asb"
        asb_generated_mode(parsed_args)
    elseif parsed_args["mode"] == "given"
        given_queries_mode(parsed_args)
    else
        throw(ArgumentError("Not implemented yet."))
    end
end