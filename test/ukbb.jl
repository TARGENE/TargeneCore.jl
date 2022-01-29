module TestsUKBB

using Test
using BGEN
using CSV
using DataFrames
using TMLEEpistasis
using CategoricalArrays
using TOML
using Serialization
using TMLE


include("helper_fns.jl")

function test_base_serialization(queryreports)
    queryreport = queryreports[1]
    @test queryreport isa TMLE.QueryReport
    test_queries((queryreport.query,),
        (Query(name="QUERY_1", case=(RSID_10 = "AG", RSID_100 = "AG"), control=(RSID_10 = "GG", RSID_100 = "GG")),)
    )
    # Check second queryreport
    queryreport = queryreports[2]
    @test queryreport isa TMLE.QueryReport
    test_queries((queryreport.query,),
        (Query(name="QUERY_2", case=(RSID_10 = "AG", RSID_100 = "AA"), control=(RSID_10 = "GG", RSID_100 = "GG")),)
    )
end


@testset "Test read_bgen function" begin
    # This file as an associated .bgi
    bgen_file = BGEN.datadir("example.8bits.bgen")
    b = TMLEEpistasis.read_bgen(bgen_file)
    @test b.idx isa BGEN.Index

    # This file as no associated .bgi
    bgen_file = BGEN.datadir("example.9bits.bgen")
    b = TMLEEpistasis.read_bgen(bgen_file)
    @test b.idx === nothing
end

@testset "Test samples_genotype" begin
    probabilities = [0.3 0.2 0.9;
                     0.5 0.2 0.05;
                     0.2 0.6 0.05]
    variant_genotypes = ["AA", "AT", "TT"]

    threshold = 0.9
    genotypes = TMLEEpistasis.samples_genotype(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2] === missing
    @test genotypes[3] == "AA"

    threshold = 0.55
    genotypes = TMLEEpistasis.samples_genotype(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1]  === missing
    @test genotypes[2] == "TT"
    @test genotypes[3] == "AA"

end

@testset "Test UKBBGenotypes function" begin
    build_query_file(threshold=0.95)
    queries = TMLEEpistasis.parse_queries(queryfile)
    genotypes = TMLEEpistasis.UKBBGenotypes(queryfile, queries)
    # I only look at the first 10 rows
    # SAMPLE_ID    
    @test genotypes[1:9, "SAMPLE_ID"] == ["sample_00$i" for i in 1:9]
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == repeat(["AG"], 10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== ["AG", "AA", "AG", missing, "AG", "AG", missing, "AG", "GG", "AG"])
    # Test column order
    @test DataFrames.names(genotypes) == ["SAMPLE_ID", "RSID_10", "RSID_100"]

    rm(queryfile)
end



@testset "Test variant_genotypes" begin
    b = TMLEEpistasis.read_bgen(BGEN.datadir("example.8bits.bgen"))
    queries = [
        Query(case=(RSID_10="AA", RSID_100="AA"), control=(RSID_10="GG", RSID_100="GG")),
        Query(case=(RSID_10="AA", RSID_100="AA"), control=(RSID_10="GG", RSID_100="GA"))
    ]
    # Defaults to "AG"
    v = variant_by_rsid(b, "RSID_10")
    @test TMLEEpistasis.variant_genotypes(v, queries) == ["AA", "AG", "GG"]
    # "GA" is specified by user so returned
    v = variant_by_rsid(b, "RSID_100")
    @test TMLEEpistasis.variant_genotypes(v, queries) == ["AA", "GA", "GG"]

end


@testset "Test preprocess" begin
    confounders = CSV.File(confoundersfile) |> DataFrame
    n = size(confounders)[1]
    genotypes = DataFrame(
        SAMPLE_ID = confounders.SAMPLE_ID, 
        RSID_1    = vcat(Vector{Missing}(missing, 10), collect(1:n-10)),
        RSID_2    = vcat(collect(1:n-10), Vector{Missing}(missing, 10))
        )

    # Continuous phenotypes
    phenotypes = CSV.File(phenotypefile, select=["eid", "continuous_phenotype"]) |> DataFrame
    T, W, y = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes;verbosity=0)

    @test y == phenotypes.continuous_phenotype[11:n-10]

    @test size(T) == (480, 2)
    @test size(W) == (480, 2)
    @test size(y) == (480,)

    @test DataFrames.names(T) == ["RSID_1", "RSID_2"]
    @test DataFrames.names(W) == ["PC1", "PC2"]

    @test T.RSID_1 == collect(1:n-20)
    @test T.RSID_1 isa CategoricalArray
    @test T.RSID_2 == collect(11:n-10)
    @test T.RSID_2 isa CategoricalArray

    # Binary phenotypes
    phenotypes = CSV.File(phenotypefile, select=["eid", "categorical_phenotype"]) |> DataFrame
    T, W, y = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes;
                verbosity=0)

    @test y isa CategoricalArray{Bool,1,UInt32}
end


@testset "Test PhenotypeTMLEEpistasis with continuous target" begin
    # Only one continuous phenotype
    estimatorfile = joinpath("config", "tmle_config.toml")
    build_query_file()
    parsed_args = Dict(
        "phenotypes" => phenotypefile,
        "confounders" => confoundersfile,
        "queries" => queryfile,
        "estimator" => estimatorfile,
        "output" => "cont_results.bin",
        "phenotypes-list" => phenotypelist_file_1,
        "verbosity" => 0,
        "adaptive-cv" => false,
        "savefull" => false
    )

    UKBBVariantRun(parsed_args)

    file = open(parsed_args["output"])

    phenotype, queryreports = deserialize(file)
    @test phenotype == :continuous_phenotype
    @test length(queryreports) == 2
    # Check queryreports
    test_base_serialization(queryreports)
    # Only the continuous phenotype has been processed
    @test eof(file)

    close(file)

    # Clean
    rm(parsed_args["output"])
    rm(parsed_args["queries"])
end


@testset "Test PhenotypeTMLEEpistasis with all phenotypes" begin
    estimatorfile = joinpath("config", "tmle_config.toml")
    build_query_file()
    parsed_args = Dict(
        "phenotypes" => phenotypefile,
        "confounders" => confoundersfile,
        "queries" => queryfile,
        "estimator" => estimatorfile,
        "output" => "full_results.csv",
        "phenotypes-list" => nothing,
        "verbosity" => 0,
        "adaptive-cv" => true,
        "savefull" => false
    )

    UKBBVariantRun(parsed_args)
    
    file = open(parsed_args["output"])
    # First phenotype
    phenotype, queryreports = deserialize(file)
    @test phenotype == :categorical_phenotype
    @test length(queryreports) == 2
    # Check queryreports
    test_base_serialization(queryreports)

    # Second phenotype
    phenotype, queryreports = deserialize(file)
    @test phenotype == :continuous_phenotype
    test_base_serialization(queryreports)

    # Only the continuous phenotype has been processed
    @test eof(file)

    close(file)

    # Clean
    rm(parsed_args["output"])
    rm(parsed_args["queries"])
end



end;

true