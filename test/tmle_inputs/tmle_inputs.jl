module TestTMLEInputs

using Test
using CSV
using DataFrames
using TargeneCore
using YAML
using StableRNGs
using BGEN
using TMLE

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
    bgen_dir = joinpath("data", "ukbb", "imputed" , "ukbb")
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
    @test_throws TargeneCore.NotAllVariantsFoundError([], variants) TargeneCore.call_genotypes(bgen_dir, variants, 0.95;)
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
    Ψ = CM(target=:toto, treatment=(A=1,), confounders=[])
    @test TargeneCore.setting_iterator(Ψ) == ((A = 1,),)
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.4) == true
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.6) == false

    Ψ = ATE(target=:toto, treatment=(A=(case=1, control=0),), confounders=[])
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

    Ψ = CM(target=:toto, treatment=(B = "CC", A = 1), confounders=[])
    @test TargeneCore.setting_iterator(Ψ) == ((A = 1, B = "CC"),)
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.1) == true
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.15) == false
    
    Ψ = ATE(target=:toto, treatment=(B=(case="AA", control="AC"), A=(case=1, control=1),), confounders=[])
    @test collect(TargeneCore.setting_iterator(Ψ)) == [(A = 1, B = "AA"), (A = 1, B = "AC")]
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.1) == true
    @test TargeneCore.satisfies_positivity(Ψ, freqs, positivity_constraint=0.2) == false
    
    Ψ = IATE(target=:toto, treatment=(B=(case="AC", control="AA"), A=(case=1, control=0),), confounders=[])
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

    Ψ = IATE(target=:toto, treatment=(B=(case="AC", control="AA"), A=(case=1, control=0), C=(control=0, case=2)), confounders=[])
    expected_settings = [
        (A = 1, B = "AC", C = 0),
        (A = 0, B = "AC", C = 0),
        (A = 1, B = "AA", C = 0),
        (A = 0, B = "AA", C = 0),
        (A = 1, B = "AC", C = 2),
        (A = 0, B = "AC", C = 2),
        (A = 1, B = "AA", C = 2),
        (A = 0, B = "AA", C = 2)]
    for (index, s) in enumerate(TargeneCore.setting_iterator(Ψ))
        @test s == expected_settings[index]
    end
end

end

true