module TestTMLEInputs

using Test
using CSV
using DataFrames
using TargeneCore
using YAML
using StableRNGs
using BGEN
using TMLE

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

#####################################################################
##################           UNIT TESTS            ##################
#####################################################################

@testset "Test genotypes_encoding" begin
    b = Bgen(BGEN.datadir("example.8bits.bgen"))
    v = variant_by_rsid(b, "RSID_10")
    minor_allele_dosage!(b, v)
    # The minor allele is the first one
    @test minor_allele(v) == alleles(v)[1]
    @test TargeneCore.genotypes_encoding(v) == ["AA", "AG", "GG"]
    # The minor allele is the second one
    v = variant_by_rsid(b, "RSID_102")
    minor_allele_dosage!(b, v)
    @test minor_allele(v) == alleles(v)[2]
    @test TargeneCore.genotypes_encoding(v) == ["AA", "AG", "GG"]
end

@testset "Test call_genotypes for a single SNP" begin
    probabilities = [-Inf 0.3 0.2 0.9 -Inf;
                    -Inf 0.5 0.2 0.05 0.3;
                    -Inf 0.2 0.6 0.05 0.7]
    variant_genotypes = [2, 1, 0]

    # No threshold results in max probability genotype being called
    threshold = nothing
    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        variant_genotypes, 
        threshold
    )
    @test all(x === y for (x,y) in zip(genotypes,[missing, 1, 0, 2, 0]))

    # If a threshold is set, the genotype is only called if the max probability is greater than the threshold
    threshold = 0.9
    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2] === genotypes[3] === missing
    @test genotypes[4] == 2
    @test genotypes[5] === missing
    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        ["AA", "AG", "GG"], 
        threshold)
    @test genotypes[4] == "AA"

    threshold = 0.55
    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2]  === missing
    @test genotypes[3] == 0
    @test genotypes[4] == 2
    @test genotypes[5] == 0

    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        ["AA", "AG", "GG"], 
        threshold)
    @test genotypes[3] == "GG"
    @test genotypes[4] == "AA"
end

@testset "Test call_genotypes for all SNPs" begin
    bgen_dir = joinpath(TESTDIR, "data", "ukbb", "imputed" , "ukbb")
    variants = Set(["RSID_10", "RSID_100", "RSID_2"])
    # With a threshold
    genotypes = TargeneCore.call_genotypes(bgen_dir, variants, 0.95)
    # I only looked at the first 10 rows
    # SAMPLE_ID    
    @test genotypes[1:9, "SAMPLE_ID"] == ["sample_00$i" for i in 1:9]
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == repeat(["AG"], 10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== ["AG", "AA", "AG", missing, "AG", "AG", missing, "AG", "GG", "AG"])
    # RSID_2: contains NaN
    @test all(genotypes[1:10, "RSID_2"] .=== [missing, "GG", missing, missing, missing, missing, missing, "GG", missing, missing])
    # Test column order
    @test DataFrames.names(genotypes) == ["SAMPLE_ID", "RSID_10", "RSID_100", "RSID_2"]

    # With no threshold
    genotypes = TargeneCore.call_genotypes(bgen_dir, variants, nothing)
    # I only looked at the first 10 rows
    # SAMPLE_ID    
    @test genotypes[1:9, "SAMPLE_ID"] == ["sample_00$i" for i in 1:9]
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == repeat(["AG"], 10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== ["AG", "AA", "AG", "AA", "AG", "AG", "AG", "AG", "GG", "AG"])
    # RSID_2: contains NaN
    @test all(genotypes[1:10, "RSID_2"] .=== [missing, "GG", "GG", "AG", "GG", "GG", "AG", "GG", "GG", "GG"])
    # Test column order
    @test DataFrames.names(genotypes) == ["SAMPLE_ID", "RSID_10", "RSID_100", "RSID_2"]

    # With missing variants -> throw ArgumentError
    variants = Set(["TOTO"])
    @test_throws TargeneCore.NotAllVariantsFoundError(variants) TargeneCore.call_genotypes(bgen_dir, variants, 0.95;)
end

end

true