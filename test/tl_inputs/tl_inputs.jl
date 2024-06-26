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
    probabilities = [NaN 0.3 0.2 0.9;
                     NaN 0.5 0.2 0.05;
                     NaN 0.2 0.6 0.05]
    variant_genotypes = [2, 1, 0]

    threshold = 0.9
    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2] === genotypes[3] === missing
    @test genotypes[4] == 2
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
    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        ["AA", "AG", "GG"], 
        threshold)
    @test genotypes[3] == "GG"
    @test genotypes[4] == "AA"
end

@testset "Test call_genotypes for all SNPs" begin
    bgen_dir = joinpath(TESTDIR, "data", "ukbb", "imputed" , "ukbb")
    variants = Set(["RSID_10", "RSID_100"])
    genotypes = TargeneCore.call_genotypes(bgen_dir, variants, 0.95)
    # I only look at the first 10 rows
    # SAMPLE_ID    
    @test genotypes[1:9, "SAMPLE_ID"] == ["sample_00$i" for i in 1:9]
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == repeat(["AG"], 10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== ["AG", "AA", "AG", missing, "AG", "AG", missing, "AG", "GG", "AG"])
    # Test column order
    @test DataFrames.names(genotypes) == ["SAMPLE_ID", "RSID_10", "RSID_100"]

    # With missing variants -> throw ArgumentError
    variants = Set(["TOTO"])
    @test_throws TargeneCore.NotAllVariantsFoundError(variants) TargeneCore.call_genotypes(bgen_dir, variants, 0.95;)
end

end

true