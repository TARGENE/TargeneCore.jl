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


@testset "Test positivity_constraint" begin
    data = DataFrame(
        A = [1, 1, 0, 1, 0, 2, 2, 1],
        B = ["AC", "CC", "AA", "AA", "AA", "AA", "AA", "AA"]
    ) 
    ## One variable
    freqs = TargeneCore.frequency_table(data, [:A])
    @test freqs == Dict(
        (A = 0,) => 0.25,
        (A = 2,) => 0.25,
        (A = 1,) => 0.5
    )
    Ψ = CM(
        outcome = :toto, 
        treatment_values = (A=1,), 
        treatment_confounders = (A=[],)
    )
    @test TargeneCore.setting_iterator(Ψ) == ((A = 1,),)
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.4) == true
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.6) == false

    Ψ = ATE(
        outcome = :toto, 
        treatment_values= (A= (case=1, control=0),), 
        treatment_confounders = (A=[],)
    )
    @test collect(TargeneCore.setting_iterator(Ψ)) == [(A = 1,), (A = 0,)]
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.2) == true
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.3) == false

    ## Two variables
    # Treatments are sorted: [:B, :A] -> [:A, :B]
    freqs = TargeneCore.frequency_table(data, [:B, :A])
    @test freqs == Dict(
        (A = 1, B = "CC") => 0.125,
        (A = 1, B = "AA") => 0.25,
        (A = 0, B = "AA") => 0.25,
        (A = 1, B = "AC") => 0.125,
        (A = 2, B = "AA") => 0.25
    )

    Ψ = CM(
        outcome = :toto, 
        treatment_values = (B = "CC", A = 1), 
        treatment_confounders = (B = [], A = [])
    )
    @test TargeneCore.setting_iterator(Ψ) == ((A = 1, B = "CC"),)
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.1) == true
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.15) == false
    
    Ψ = ATE(
        outcome = :toto, 
        treatment_values = (B=(case="AA", control="AC"), A=(case=1, control=1),), 
        treatment_confounders = (B = (), A = (),)
    )
    @test collect(TargeneCore.setting_iterator(Ψ)) == [(A = 1, B = "AA"), (A = 1, B = "AC")]
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.1) == true
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.2) == false
    
    Ψ = IATE(
        outcome = :toto, 
        treatment_values = (B=(case="AC", control="AA"), A=(case=1, control=0),), 
        treatment_confounders = (B=(), A=()), 
    )
    @test collect(TargeneCore.setting_iterator(Ψ)) == [
        (A = 1, B = "AC")  (A = 1, B = "AA")
        (A = 0, B = "AC")  (A = 0, B = "AA")]
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=1.) == false
    freqs = Dict(
        (A = 1, B = "CC") => 0.125,
        (A = 1, B = "AA") => 0.25,
        (A = 0, B = "AA") => 0.25,
        (A = 0, B = "AC") => 0.25,
        (A = 1, B = "AC") => 0.125,
        (A = 2, B = "AA") => 0.25
    )
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.3) == false
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.1) == true

    Ψ = IATE(
        outcome = :toto, 
        treatment_values = (B=(case="AC", control="AA"), A=(case=1, control=0), C=(control=0, case=2)), 
        treatment_confounders = (B=(), A=(), C=())
    )
    expected_settings = Set([
        (A = 1, B = "AC", C = 0),
        (A = 0, B = "AC", C = 0),
        (A = 1, B = "AA", C = 0),
        (A = 0, B = "AA", C = 0),
        (A = 1, B = "AC", C = 2),
        (A = 0, B = "AC", C = 2),
        (A = 1, B = "AA", C = 2),
        (A = 0, B = "AA", C = 2)])
    @test expected_settings == Set(TargeneCore.setting_iterator(Ψ))
end

end

true