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


@testset "Test UKBBGenotypes function" begin
    build_query_file(threshold=0.95)
    genotypes = GenesInteraction.UKBBGenotypes(queryfile)
    # I only looked at the firs 10 rows
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == repeat(["AG"], 10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== ["AG", "AA", "AG", missing, "AG", "AG", missing, "AG", "GG", "AG"])

    rm(queryfile)
end


end;

true