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
    genotypes = GenesInteraction.UKBBGenotypes(queryfile, queries["QUERY_1"])
    # I only looked at the firs 10 rows
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
        "QUERY_1" => (RSID_10 = ["AG", "GG"], RSID_100 = ["AA", "GG"]),
        "QUERY_2" => (RSID_10 = ["AA", "GG"], RSID_100 = ["AA", "GG"])
    )
    rm(queryfile)
end


@testset "Test variant_genotypes" begin
    b = GenesInteraction.read_bgen(BGEN.datadir("example.8bits.bgen"))
    v = variant_by_rsid(b, "RSID_10")

    query = (RSID_10=["AG", "GG"],)
    @test GenesInteraction.variant_genotypes(v, query) == ["AA", "AG", "GG"]

    query = (RSID_10=["GA", "GG"],)
    @test GenesInteraction.variant_genotypes(v, query) == ["AA", "GA", "GG"]

    query = (RSID_10=["AA", "GG"],)
    @test GenesInteraction.variant_genotypes(v, query) == ["AA", "AG", "GG"]
end


# @testset "Epistasis estimation run" begin
#     tmle_config = joinpath("config", "tmle_continuous.toml")
#     build_query_file()
#     outfile = "tmleepistasis_results.csv"
#     TMLEEpistasisUKBB(genotypefile, 
#                     phenotypefile, 
#                     confoundersfile, 
#                     snpfile,
#                     tmle_config,
#                     outfile)
    
#     try
#         results = CSV.File(outfile) |> DataFrame
#         @test results[!, "Estimate"] isa Vector
#         @test results[!, "Standard Error"] isa Vector
#     finally
#         rm(outfile)
#         rm(queryfile)
#     end

# end

end;

true