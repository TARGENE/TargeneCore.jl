module TestUtils

using Test
using TargeneCore
using BGEN
using DataFrames

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "testutils.jl"))

@testset "Test Estimands Accessors" begin
    traits = TargeneCore.read_csv_file(joinpath(TESTDIR, "data", "traits_1.csv"))
    pcs = TargeneCore.load_flash_pca_results(joinpath(TESTDIR, "data", "pcs.csv"))
    # extraW, extraT, extraC are parsed from all estimands_files
    estimands = make_estimands_configuration().estimands
    # get_treatments, get_outcome, ...
    ## Simple Estimand
    Ψ = estimands[1]
    @test get_outcome(Ψ) == :ALL
    @test get_treatments(Ψ) == (:RSID_2, :TREAT_1)
    @test get_all_confounders(Ψ) == ()
    @test get_outcome_extra_covariates(Ψ) == ()
    ## Simple Estimand with Confounders
    Ψ = estimands[2]
    @test get_outcome(Ψ) == :ALL
    @test get_treatments(Ψ) == (:RSID_2,)
    @test get_all_confounders(Ψ) == (Symbol("22001"),)
    @test get_confounders(Ψ, :RSID_2) == (Symbol("22001"),)
    @test get_outcome_extra_covariates(Ψ) == (Symbol("21003"), :COV_1)
    ## JointEstimand
    Ψ = estimands[5]
    @test get_outcome(Ψ) == :ALL
    @test get_treatments(Ψ) == (:RSID_198, :RSID_2)
    @test get_all_confounders(Ψ) == (:PC1, :PC2)
    @test get_confounders(Ψ, :RSID_198) == (:PC2,)
    @test get_confounders(Ψ, :RSID_2) == (:PC1, )
    @test get_outcome_extra_covariates(Ψ) == (Symbol("22001"), )
    ## Bad JointEstimand
    Ψ = JointEstimand(
        CM(
            outcome = "Y1",
            treatment_values = (RSID_3 = "GG", RSID_198 = "AG"),
            treatment_confounders = (RSID_3 = [], RSID_198 = []),
            outcome_extra_covariates = [22001]
        ),
        CM(
            outcome = "Y2",
            treatment_values = (RSID_2 = "AA", RSID_198 = "AG"),
            treatment_confounders = (RSID_2 = [:PC1], RSID_198 = []),
            outcome_extra_covariates = []
        )
    )
    @test_throws ArgumentError get_outcome(Ψ)
    @test_throws ArgumentError get_treatments(Ψ)
    @test_throws ArgumentError get_all_confounders(Ψ)
    @test_throws ArgumentError get_outcome_extra_covariates(Ψ)
    # get_variables
    variables = TargeneCore.get_variables(estimands, traits, pcs)
    @test variables.genetic_variants == Set([:RSID_198, :RSID_2])
    @test variables.outcomes == Set([:BINARY_1, :CONTINUOUS_2, :CONTINUOUS_1, :BINARY_2])
    @test variables.pcs == Set([:PC1, :PC2])
end

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

@testset "Test pvalue_or_nan" begin
    estimates = make_estimates()

    @test TargeneCore.pvalue_or_nan(estimates[1].TMLE) isa Float64
    @test TargeneCore.pvalue_or_nan(estimates[3].TMLE) isa Float64
    @test TargeneCore.pvalue_or_nan(estimates[5].TMLE) === NaN

    @test TargeneCore.pvalue_or_nan(estimates[1].TMLE) ≠ TargeneCore.pvalue_or_nan(estimates[1].TMLE, -1)
    @test TargeneCore.pvalue_or_nan(estimates[3].TMLE) ≠ TargeneCore.pvalue_or_nan(estimates[3].TMLE, [-1., -1.])
end

end

true