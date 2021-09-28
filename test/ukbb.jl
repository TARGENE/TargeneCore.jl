module TestsUKBB

using Test
using BGEN
using CSV
using DataFrames
using GenesInteraction
using CategoricalArrays

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

@testset "Test variant_alleles_indices function" begin
    bgen_file = BGEN.datadir("example.8bits.bgen")
    b = GenesInteraction.read_bgen(bgen_file)
    v = variant_by_rsid(b, "RSID_10")
    # A prior call to probabilities! is required to load the data
    # in order to see if it is phased or not
    probabilities!(b, v)

    @test GenesInteraction.variant_alleles_indices(v, "A/A") == 1
    @test GenesInteraction.variant_alleles_indices(v, "A/G") == 2
    @test GenesInteraction.variant_alleles_indices(v, "G/A") == 2
    @test GenesInteraction.variant_alleles_indices(v, "G/G") == 3
end


@testset "Test retrieve_alleles function" begin
    probabilities = [0.2 0.8 0.1;
                     0.7 0.1 0.1;
                     0.1 0.1 0.8]

    @test all(GenesInteraction.retrieve_alleles(probabilities, 2, 3; threshold=0.79) .=== [missing, missing, false])
    @test all(GenesInteraction.retrieve_alleles(probabilities, 2, 3; threshold=0.69) .=== [true, missing, false])
    @test all(GenesInteraction.retrieve_alleles(probabilities, 2, 1; threshold=0.69) .=== [true, false, missing])
    @test all(GenesInteraction.retrieve_alleles(probabilities, 2, 1; threshold=0.79) .=== [missing, false, missing])
    @test all(GenesInteraction.retrieve_alleles(probabilities, 3, 1; threshold=0.69) .=== [missing, false, true])
    @test all(GenesInteraction.retrieve_alleles(probabilities, 3, 1; threshold=0.79) .=== [missing, false, true])
end


@testset "Test UKBBGenotypes function" begin
    queryfile = joinpath("data", "query.csv")
    build_query_file(queryfile)
    
    genotypes = GenesInteraction.UKBBGenotypes(queryfile; threshold=0.95)
    #I only looked at the firs 10 rows
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == [true, true, true, true, true, true, true, true, true, true]
    # RSID_100
    # More varied 
    @test all(genotypes[1:10, "RSID_100"] .=== [missing, false, missing, missing, missing, missing, missing, missing, true, missing])

    rm(queryfile)
end


end;

true