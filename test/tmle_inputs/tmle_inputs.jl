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
    @test TargeneCore.genotypes_encoding(v) == [2, 1, 0]
    @test TargeneCore.genotypes_encoding(v, asint=false) == ["AA", "AG", "GG"]
    # The minor allele is the second one
    v = variant_by_rsid(b, "RSID_102")
    minor_allele_dosage!(b, v)
    @test minor_allele(v) == alleles(v)[2]
    @test TargeneCore.genotypes_encoding(v) == [0, 1, 2]
    @test TargeneCore.genotypes_encoding(v, asint=false) == ["AA", "AG", "GG"]
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
    snp_list = ["RSID_10", "RSID_100"]
    ## With asint=true
    genotypes = TargeneCore.call_genotypes(bgen_dir, snp_list, 0.95; asint=true)
    # I only look at the first 10 rows
    # SAMPLE_ID    
    @test genotypes[1:9, "SAMPLE_ID"] == ["sample_00$i" for i in 1:9]
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == ones(10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== [1, 2, 1, missing, 1, 1, missing, 1, 0, 1])
    # Test column order
    @test DataFrames.names(genotypes) == ["SAMPLE_ID", "RSID_10", "RSID_100"]

    ## With asint=false
    genotypes = TargeneCore.call_genotypes(bgen_dir, snp_list, 0.95; asint=false)
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
    snp_list = ["TOTO"]
    @test_throws TargeneCore.NotAllVariantsFoundError([], snp_list) TargeneCore.call_genotypes(bgen_dir, snp_list, 0.95;)
end

@testset "Test treatment settings generation & positivity_respecting_parameters" begin
    freqs = Dict(
        ("A", 1) => 0.01,
        ("B", 1) => 0.02,
        ("A", 0) => 0.03,
        ("B", 0) => 0.04,
        ("C", 1) => 0.05,
        ("C", 0) => 0.06,
    )
    ## Multiple variables: IATE, ATE. CM
    # IATE
    # All 4 product terms A/0, A/1, B/0, B/1 are checked
    interaction_setting = (["A", "B"], [0, 1])
    @test collect(TargeneCore.setting_iterator(IATE, interaction_setting)) == 
        [("A", 0)  ("A", 1)
        ("B", 0)  ("B", 1)]
    @test TargeneCore.satisfies_positivity(IATE, interaction_setting, freqs; positivity_constraint=0.01) == true
    interaction_setting = (["B", "C"], [0, 1])
    @test TargeneCore.satisfies_positivity(IATE, interaction_setting, freqs; positivity_constraint=0.02) == true
    @test TargeneCore.satisfies_positivity(IATE, interaction_setting, freqs; positivity_constraint=0.03) == false

    # ATE
    # Only look at A/0 and B/1
    ate_setting = (["A", "B"], [0, 1])
    @test collect(TargeneCore.setting_iterator(ATE, ate_setting)) == 
        [("A", 0), ("B", 1)]
    @test TargeneCore.satisfies_positivity(ATE, ate_setting, freqs; positivity_constraint=0.03) == false
    @test TargeneCore.satisfies_positivity(ATE, ate_setting, freqs; positivity_constraint=0.01) == true

    # CM
    cm_setting = (["A",], [0,])
    @test collect(TargeneCore.setting_iterator(CM, cm_setting)) == [("A", 0)]
    @test TargeneCore.satisfies_positivity(CM, cm_setting, freqs; positivity_constraint=0.04) == false
    @test TargeneCore.satisfies_positivity(CM, cm_setting, freqs; positivity_constraint=0.03) == true


    ## One variable: ATE, CM
    # Only look at A/0 and B/1
    freqs = Dict(
        ("A",) => 0.01,
        ("B",) => 0.02,
    )
    # ATE
    ate_setting = (["A", "B"],)
    @test collect(TargeneCore.setting_iterator(ATE, ate_setting)) == 
        [("A",), ("B", )]
    @test TargeneCore.satisfies_positivity(ATE, ate_setting, freqs; positivity_constraint=0.03) == false
    @test TargeneCore.satisfies_positivity(ATE, ate_setting, freqs; positivity_constraint=0.01) == true
    # CM
    cm_setting = (["A"],)
    @test TargeneCore.satisfies_positivity(ATE, ate_setting, freqs; positivity_constraint=0.02) == false
    @test TargeneCore.satisfies_positivity(ATE, ate_setting, freqs; positivity_constraint=0.01) == true
end

@testset "Test frequency_table" begin
    data = DataFrame(
        A = [1, 1, 0, 1, 0, 2, 2, 1],
        B = ["AC", "CC", "AA", "AA", "AA", "AA", "AA", "AA"]
    ) 
    freqs = TargeneCore.frequency_table(data, [:A])
    @test freqs == Dict(
        (0,) => 0.25,
        (2,) => 0.25,
        (1,) => 0.5
    )
    freqs = TargeneCore.frequency_table(data, [:A, :B])
    @test freqs == Dict(
        (1, "CC") => 0.125,
        (1, "AA") => 0.25,
        (0, "AA") => 0.25,
        (1, "AC") => 0.125,
        (2, "AA") => 0.25
    )

end

end

true