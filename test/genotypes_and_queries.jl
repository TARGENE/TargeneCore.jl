module TestGenerateQueries

using Test
using TMLEEpistasis
using DataFrames
using CSV
using TOML
using BGEN


function write_queries()
    open("test_query_1.toml", "w") do io
        TOML.print(io, Dict(
            "SNPS" => Dict("RSID_2" => "chr12"),
            "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR" => Dict("RSID_2"=>"GG -> GA")
        ))
    end
    open("test_query_2.toml", "w") do io
        TOML.print(io, Dict(
            "SNPS" => Dict("RSID_2" => "chr12", "RSID_198" => "chr12"),
            "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR" => Dict("RSID_2"=>"GG -> GA", "RSID_198"=>"GG -> GA")
        ))
    end
end

function remove_queries()
    rm("test_query_1.toml")
    rm("test_query_2.toml")
end

@testset "Test allele_specific_binding_snps" begin
    bQTLs = TMLEEpistasis.allele_specific_binding_snps(joinpath("data", "asb_files", "asb_"))
    @test bQTLs == DataFrame(
        ID=["RSID_17", "RSID_99", "RSID_198"],
        CHROM=["chr12", "chr12", "chr12"],
        QTL_TYPE=[:bQTL, :bQTL, :bQTL],
    )
end

@testset "Test trans_actors(path)" begin
    eQTLs = TMLEEpistasis.trans_actors(joinpath("data", "trans_actors_fake.csv"))
    @test eQTLs == DataFrame(
        ID=["RSID_102", "RSID_2"], 
        CHROM=["chr12", "chr12"], 
        QTL_TYPE=[:eQTL, :eQTL]
    )
end

@testset "Test check_chr_format" begin
    TMLEEpistasis.check_chr_format(["chr1", "chr43"])
    @test_throws AssertionError TMLEEpistasis.check_chr_format(["chr", "chr43"])
    @test_throws AssertionError TMLEEpistasis.check_chr_format(["12", "chr43"])
end

@testset "Test exclude" begin
    snps = DataFrame(
        ID=["RSID_17", "RSID_12"], 
        CHROM=["chr12", "chr13"]
        )

    TMLEEpistasis.exclude!(snps, nothing)
    @test size(snps, 1) == 2

    TMLEEpistasis.exclude!(snps, joinpath("data", "pb_snps.csv"))
    @test size(snps, 1) == 1
end

@testset "Test chr_filenames" begin
    # Missing chromosome
    chr_prefix = joinpath("data", "ukbb", "imputed" ,"ukbb")
    snps = DataFrame(
        ID=["RSID_17", "RSID_12"], 
        CHROM=["chr12", "chr13"]
        )
    @test_throws AssertionError TMLEEpistasis.add_chrfiles(snps, chr_prefix)
    # Everything fine
    snps = DataFrame(
        ID=["RSID_17", "RSID_12"], 
        CHROM=["chr12", "chr12"]
        )
    snps = TMLEEpistasis.add_chrfiles(snps, chr_prefix)
    @test snps.ID == ["RSID_17", "RSID_12"]
    @test all(x isa Bgen for x in snps.BGEN_FILE)
end

@testset "Test variant_genotypes" begin
    b = Bgen(BGEN.datadir("example.8bits.bgen"))
    v = variant_by_rsid(b, "RSID_10")
    minor_allele_dosage!(b, v)
    major = major_allele(v)
    minor = minor_allele(v)
    @test TMLEEpistasis.genotypes(v, minor, major) == ["AA", "GA", "GG"]
end

@testset "Test impute_genotypes for a single SNP" begin
    probabilities = [0.3 0.2 0.9;
                     0.5 0.2 0.05;
                     0.2 0.6 0.05]
    variant_genotypes = ["AA", "AT", "TT"]

    threshold = 0.9
    genotypes = TMLEEpistasis.impute_genotypes(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2] === missing
    @test genotypes[3] == "AA"

    threshold = 0.55
    genotypes = TMLEEpistasis.impute_genotypes(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1]  === missing
    @test genotypes[2] == "TT"
    @test genotypes[3] == "AA"
end

@testset "Test impute_genotypes for all SNPs" begin
    bgenfile = Bgen(joinpath("data", "ukbb", "imputed" , "ukbb_chr12.bgen"))
    snps = DataFrame(
        ID=["RSID_10", "RSID_100"], 
        CHROM=["chr12", "chr12"], 
        BGEN_FILE=[bgenfile, bgenfile]
        )
    genotypes, rsid_to_minor_major = TMLEEpistasis.impute_genotypes(snps, 0.95)
    # I only look at the first 10 rows
    # SAMPLE_ID    
    @test genotypes[1:9, "SAMPLE_ID"] == ["sample_00$i" for i in 1:9]
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == repeat(["GA"], 10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== ["GA", "AA", "GA", missing, "GA", "GA", missing, "GA", "GG", "GA"])
    # Test column order
    @test DataFrames.names(genotypes) == ["SAMPLE_ID", "RSID_10", "RSID_100"]
    # Check rsid_to_minor_major
    @test rsid_to_minor_major == Dict(
        "RSID_10"  => (minor = "A", major = "G"),
        "RSID_100" => (minor = "A", major = "G")
    )
end

@testset "Test build_queries" begin
    eQTLs = DataFrame(ID=["RSID_1", "RSID_2"], CHROM=["chr12", "chr1"])
    bQTLs = DataFrame(ID=["RSID_3"], CHROM=["chr4"])
    genotypes = DataFrame(
        RSID_1=["AA", "AG", "GG", "GG", "AG", "AA"],
        RSID_2=["CC", "CG", "GG", "CC", "GG", "CC"],
        RSID_3=["CC", "GG", "GG", "GC", "GG", "CC"]
        )
    rsid_to_major_minor = Dict(
        "RSID_1" => (minor="G", major="A"),
        "RSID_2" => (minor="G", major="C"),
        "RSID_3" => (minor="C", major="G"),
    )
    
    # no occurence of some genotypes interactions
    @test TMLEEpistasis.build_queries(eQTLs, bQTLs, genotypes, rsid_to_major_minor, 0.) == []

    # one interaction ok min freq is 1/7
    genotypes = DataFrame(
        RSID_1=["AA", "AA", "GG", "GG", "AG", "AA", "AG"],
        RSID_2=["CC", "CG", "GG", "CC", "GG", "CC", "CC"],
        RSID_3=["CC", "GG", "GG", "GC", "GC", "GC", "GG"]
        )
    queries = TMLEEpistasis.build_queries(eQTLs, bQTLs, genotypes, rsid_to_major_minor, 0.1)
    query = only(queries)
    @test query == Dict(
        "SNPS" => Dict("RSID_1" => "chr12", "RSID_3" => "chr4"),
        "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR" => Dict("RSID_1"=>"AA -> AG", "RSID_3"=>"GG -> GC")
        )
    # chaning the threshold removes the interaction
    @test TMLEEpistasis.build_queries(eQTLs, bQTLs, genotypes, rsid_to_major_minor, 0.2) == []
end

@testset "Test read_queries and snps_from_queries" begin
    write_queries()
    queries = TMLEEpistasis.read_queries("test_query")
    @test size(queries, 1) == 2
    @test all(haskey(q, "SNPS") for q in queries)
    @test all(haskey(q, "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR") for q in queries)

    snps = TMLEEpistasis.snps_from_queries(queries)

    @test snps == DataFrame(
        ID=["RSID_2", "RSID_198"],
        CHROM=["chr12", "chr12"]
    )
    remove_queries()
end

@testset "Test build_genotypes_and_queries: asb_generated_mode" begin
    outdir = joinpath("data", "outdir")

    mkdir(outdir)

    parsed_args = Dict(
        "mode" => "asb",
        "asb-prefix" => joinpath("data", "asb_files", "asb_"),
        "outdir" => outdir,
        "call-threshold" => 0.8,
        "chr-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"),
        "trans-actors" => joinpath("data", "trans_actors_fake.csv"),
        "exclude" => joinpath("data", "pb_snps.csv"),
        "minor-genotype-freq" => 0.001
    )

    build_genotypes_and_queries(parsed_args)
    # 500 individuals, 1 sample_id, 2 bQTLs, 2 eQTLs
    genotypes = CSV.read(joinpath(outdir, "genotypes.csv"), DataFrame)
    @test size(genotypes) == (500, 5)
    # queries
    query₁ = TOML.parse(open(joinpath(outdir, "query_RSID_2__RSID_99.toml")))
    @test query₁ == Dict(
        "SNPS" => Dict("RSID_2"=>"chr12", "RSID_99"=>"chr12"),
        "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR" => Dict("RSID_2"=>"GG -> GA", "RSID_99"=>"GG -> GA")
    )
    query₂ = TOML.parse(open(joinpath(outdir, "query_RSID_2__RSID_198.toml")))
    @test query₂ == Dict(
        "SNPS" => Dict("RSID_2"=>"chr12", "RSID_198"=>"chr12"),
        "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR" => Dict("RSID_2"=>"GG -> GA", "RSID_198"=>"GG -> GA")
    )
    query₃ = TOML.parse(open(joinpath(outdir, "query_RSID_102__RSID_99.toml")))
    @test query₃ == Dict(
        "SNPS" => Dict("RSID_102"=>"chr12", "RSID_99"=>"chr12"),
        "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR" => Dict("RSID_102"=>"AA -> AG", "RSID_99"=>"GG -> GA")
    )
    query₄ = TOML.parse(open(joinpath(outdir, "query_RSID_102__RSID_198.toml")))
    @test query₄ == Dict(
        "SNPS" => Dict("RSID_102"=>"chr12", "RSID_198"=>"chr12"),
        "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR" => Dict("RSID_102"=>"AA -> AG", "RSID_198"=>"GG -> GA")
    )

    # Increase filter
    parsed_args["minor-genotype-freq"] = 0.1
    rm(outdir; recursive=true)
    mkdir(outdir)
    build_genotypes_and_queries(parsed_args)
    # 500 individuals, 1 sample_id, 2 bQTLs, 2 eQTLs
    genotypes = CSV.read(joinpath(outdir, "genotypes.csv"), DataFrame)
    @test size(genotypes) == (500, 5)
    # no query on this filter
    @test size(readdir(outdir), 1) == 1

    # Cleanup
    rm(outdir; recursive=true)
end


@testset "Test build_genotypes_and_queries:queries given" begin
    outdir = joinpath("data", "outdir")

    mkdir(outdir)
    write_queries()
    parsed_args = Dict(
        "mode" => "given",
        "query-prefix" => "test_query_",
        "outdir" => outdir,
        "call-threshold" => 0.8,
        "chr-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"),
        "exclude" => joinpath("data", "pb_snps.csv"),
    )

    build_genotypes_and_queries(parsed_args)

    genotypes = CSV.read(joinpath(outdir, "genotypes.csv"), DataFrame)
    @test size(genotypes) == (500, 3)
    @test names(genotypes) == ["SAMPLE_ID", "RSID_2", "RSID_198"]

    query₁ = TOML.parse(open(joinpath(outdir, "query_RSID_2.toml")))
    @test query₁ == Dict(
        "SNPS" => Dict("RSID_2"=>"chr12"),
        "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR" => Dict("RSID_2"=>"GG -> GA")
    )

    query₂ = TOML.parse(open(joinpath(outdir, "query_RSID_2__RSID_198.toml")))
    @test query₂ == Dict(
        "SNPS" => Dict("RSID_2"=>"chr12", "RSID_198"=>"chr12"),
        "HOMOZYGOUS_MAJOR_TO_MAJOR_MINOR" => Dict("RSID_2"=>"GG -> GA", "RSID_198"=>"GG -> GA")
    )

    rm(outdir; recursive=true)
    remove_queries()

end

end;

true