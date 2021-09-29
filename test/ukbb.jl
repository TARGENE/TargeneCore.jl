module TestsUKBB

using Test
using BGEN
using CSV
using DataFrames
using GenesInteraction
using CategoricalArrays
using TOML

include("helper_fns.jl")

@testset "Test read_bgen function" begin
    # This file as an associated .bgi
    bgen_file = BGEN.datadir("example.8bits.bgen")
    b = GenesInteraction.read_bgen(bgen_file)
    @test b.idx isa BGEN.Index

    # This file as no associated .bgi
    bgen_file = BGEN.datadir("example.9bits.bgen")
    b = GenesInteraction.read_bgen(bgen_file)
    @test b.idx === nothing
end

@testset "Test samples_genotype" begin
    probabilities = [0.3 0.2 0.9;
                     0.5 0.2 0.05;
                     0.2 0.6 0.05]
    variant_genotypes = ["AA", "AT", "TT"]

    threshold = 0.9
    genotypes = GenesInteraction.samples_genotype(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2] === missing
    @test genotypes[3] == "AA"

    threshold = 0.55
    genotypes = GenesInteraction.samples_genotype(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1]  === missing
    @test genotypes[2] == "TT"
    @test genotypes[3] == "AA"

end

@testset "Test UKBBGenotypes function" begin
    build_query_file(threshold=0.95)
    queries = GenesInteraction.parse_queries(queryfile)
    genotypes = GenesInteraction.UKBBGenotypes(queryfile, queries)
    # I only look at the first 10 rows
    # SAMPLE_ID    
    @test genotypes[1:9, "SAMPLE_ID"] == ["sample_00$i" for i in 1:9]
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == repeat(["AG"], 10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== ["AG", "AA", "AG", missing, "AG", "AG", missing, "AG", "GG", "AG"])

    rm(queryfile)
end

@testset "Test parse_queries" begin
    build_query_file()
    queries = GenesInteraction.parse_queries(queryfile)
    @test queries == Dict(
        "QUERY_1" => (RSID_10 = ["AG", "GG"], RSID_100 = ["AG", "GG"]),
        "QUERY_2" => (RSID_10 = ["AG", "GG"], RSID_100 = ["AA", "GG"])
    )
    rm(queryfile)
end


@testset "Test variant_genotypes" begin
    b = GenesInteraction.read_bgen(BGEN.datadir("example.8bits.bgen"))
    queries = Dict(
        "query_1" => (RSID_10=["AA", "GG"], RSID_100=["AA", "GG"]),
        "query_2" => (RSID_10=["AA", "GG"], RSID_100=["AA", "GA"])
    )
    # Defaults to "AG"
    v = variant_by_rsid(b, "RSID_10")
    @test GenesInteraction.variant_genotypes(v, queries) == ["AA", "AG", "GG"]
    # "GA" is specified by user so returned
    v = variant_by_rsid(b, "RSID_100")
    @test GenesInteraction.variant_genotypes(v, queries) == ["AA", "GA", "GG"]

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
    phenotypes = CSV.File(continuous_phenotype_file) |> DataFrame
    T, W, y = GenesInteraction.preprocess(genotypes, confounders, phenotypes;
                typeoftarget="continuous", verbosity=0)

    @test y == phenotypes.phenotype[11:n-10]

    @test size(T) == (480, 2)
    @test size(W) == (480, 2)
    @test size(y) == (480,)

    @test names(T) == ["RSID_1", "RSID_2"]
    @test names(W) == ["PC1", "PC2"]

    @test T.RSID_1 == collect(1:n-20)
    @test T.RSID_1 isa CategoricalArray
    @test T.RSID_2 == collect(11:n-10)
    @test T.RSID_2 isa CategoricalArray

    # Binary phenotypes
    phenotypes = CSV.File(categorical_phenotype_file) |> DataFrame
    T, W, y = GenesInteraction.preprocess(genotypes, confounders, phenotypes;
                typeoftarget="binary",
                verbosity=0)

    @test y isa CategoricalArray{Bool,1,UInt32}
end


@testset "Test TMLEEpistasisUKBB with continuous target" begin
    estimatorfile = joinpath("config", "tmle_continuous.toml")
    build_query_file()
    outfile = "cont_results.csv"
    TMLEEpistasisUKBB(continuous_phenotype_file, 
                      confoundersfile, 
                      queryfile, 
                      estimatorfile,
                      outfile;
                      verbosity=0)
    
    results = CSV.File(outfile) |> DataFrame
    @test results.QUERY == ["QUERY_1", "QUERY_2"]
    @test names(results) == ["QUERY", "ESTIMATE", "PVALUE", "LOWER_BOUND", "UPPER_BOUND", "STD_ERROR"]
    @test size(results) == (2, 6)
    rm(outfile)
    rm(queryfile)
end


@testset "Test TMLEEpistasisUKBB with categorical target" begin
    estimatorfile = joinpath("config", "tmle_categorical.toml")
    build_query_file()
    outfile = "cat_results.csv"
    TMLEEpistasisUKBB(categorical_phenotype_file, 
                      confoundersfile, 
                      queryfile, 
                      estimatorfile,
                      outfile;
                      verbosity=0)
    
    results = CSV.File(outfile) |> DataFrame
    @test results.QUERY == ["QUERY_1", "QUERY_2"]
    @test names(results) == ["QUERY", "ESTIMATE", "PVALUE", "LOWER_BOUND", "UPPER_BOUND", "STD_ERROR"]
    @test size(results) == (2, 6)
    rm(outfile)
    rm(queryfile)
end


end;

true