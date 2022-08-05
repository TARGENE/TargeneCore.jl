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
                DataFrames.select(filter(x -> x.isASB === true && x.CHROM ∉ ("chrX", "chrY"), new_bQTLs), [:ID, :CHROM])
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


exclude!(snps, exclusion_file::Nothing) = nothing
exclude!(snps, exclusion_file) = filter!(:ID => x -> x ∉ readlines(exclusion_file), snps)

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

            counts = DataFrames.combine(
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

function write_outputs(genotypes, queries, parsed_args)
    # Genotypes
    sample_ids = Set(CSV.read(parsed_args["sample-ids"], DataFrame, header=false, types=String)[!, 1])
    filtered_genotypes = filter(x -> x.SAMPLE_ID ∈ sample_ids, genotypes)
    CSV.write(joinpath(parsed_args["outdir"], "genotypes.csv"), filtered_genotypes)
    # Queries
    for query in queries
        query_filename = "query_"
        for (snp, _) in query["SNPS"]
            query_filename = string(query_filename, snp, "__")
        end
        query_filename = string(query_filename[1:end-2], ".toml")
        open(joinpath(parsed_args["outdir"], query_filename), "w") do io
            TOML.print(io, query)
        end
    end
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
    write_outputs(genotypes, queries, parsed_args)
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
    write_outputs(genotypes, queries, parsed_args)
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