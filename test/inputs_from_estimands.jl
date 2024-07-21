module TestFromParamFiles

using Test
using CSV
using DataFrames
using TargeneCore
using YAML
using TMLE
using Arrow
using Serialization
using Random
using StableRNGs

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "testutils.jl"))

#####################################################################
###############               UNIT TESTS              ###############
#####################################################################

@testset "Test get_variables" begin
    traits = TargeneCore.read_csv_file(joinpath(TESTDIR, "data", "traits_1.csv"))
    pcs = TargeneCore.read_csv_file(joinpath(TESTDIR, "data", "pcs.csv"))
    # extraW, extraT, extraC are parsed from all estimands_files
    estimands = make_estimands_configuration().estimands
    # get_treatments, get_outcome, ...
    ## Simple Estimand
    Ψ = estimands[1]
    @test TargeneCore.get_outcome(Ψ) == :ALL
    @test TargeneCore.get_treatments(Ψ) == keys(Ψ.treatment_values)
    @test TargeneCore.get_all_confounders(Ψ) == ()
    @test TargeneCore.get_outcome_extra_covariates(Ψ) == ()
    ## Simple Estimand with Confounders
    Ψ = estimands[2]
    @test TargeneCore.get_outcome(Ψ) == :ALL
    @test TargeneCore.get_treatments(Ψ) == keys(Ψ.treatment_values)
    @test TargeneCore.get_all_confounders(Ψ) == (Symbol("22001"),)
    @test TargeneCore.get_confounders(Ψ, :RSID_2) == (Symbol("22001"),)
    @test TargeneCore.get_outcome_extra_covariates(Ψ) == (Symbol("21003"), :COV_1)
    ## ComposedEstimand
    Ψ = estimands[5]
    @test TargeneCore.get_outcome(Ψ) == :ALL
    @test TargeneCore.get_treatments(Ψ) == keys(Ψ.args[1].treatment_values)
    @test TargeneCore.get_all_confounders(Ψ) == (:PC1, :PC2)
    @test TargeneCore.get_confounders(Ψ, :RSID_198) == (:PC2,)
    TargeneCore.get_confounders(Ψ, :RSID_2) == (:PC1, )
    @test TargeneCore.get_outcome_extra_covariates(Ψ) == (Symbol("22001"), )
    ## Bad ComposedEstimand
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
    @test_throws ArgumentError TargeneCore.get_outcome(Ψ)
    @test_throws ArgumentError TargeneCore.get_treatments(Ψ)
    @test_throws ArgumentError TargeneCore.get_all_confounders(Ψ)
    @test_throws ArgumentError TargeneCore.get_outcome_extra_covariates(Ψ)
    # get_variables
    variables = TargeneCore.get_variables(estimands, traits, pcs)
    @test variables.genetic_variants == Set([:RSID_198, :RSID_2])
    @test variables.outcomes == Set([:BINARY_1, :CONTINUOUS_2, :CONTINUOUS_1, :BINARY_2])
    @test variables.pcs == Set([:PC1, :PC2])
end

@testset "Test adjust_estimand_fields" begin
    genotypes = DataFrame(
        SAMPLE_ID = [1, 2, 3],
        RSID_198 = ["GA", "GG", "AA"],
        RSID_2 = ["AG", "GG", "AA"],
    )
    pcs = Set([:PC1, :PC2])
    variants_alleles = Dict(:RSID_198 => Set(genotypes.RSID_198))
    estimands = make_estimands_configuration().estimands
    # RS198 AG is not in the genotypes but GA is
    Ψ = estimands[4]
    @test Ψ.treatment_values.RSID_198 == (case="AG", control="AA")
    new_Ψ = TargeneCore.adjust_estimand_fields(Ψ, variants_alleles, pcs)
    @test new_Ψ.outcome == Ψ.outcome
    @test new_Ψ.outcome_extra_covariates == Ψ.outcome_extra_covariates
    @test new_Ψ.treatment_confounders == (RSID_198 = (:PC1, :PC2), RSID_2 = (:PC1, :PC2))
    @test new_Ψ.treatment_values == (
        RSID_198 = (case = "GA", control = "AA"),
        RSID_2 = (case = "AA", control = "GG")
    )

    # JointEstimand
    Ψ = estimands[5]
    @test Ψ.args[1].treatment_values == (RSID_198 = "AG", RSID_2 = "GG")
    @test Ψ.args[2].treatment_values == (RSID_198 = "AG", RSID_2 = "AA")
    new_Ψ = TargeneCore.adjust_estimand_fields(Ψ, variants_alleles, pcs)
    for index in 1:length(Ψ.args)
        @test new_Ψ.args[index].outcome == Ψ.args[index].outcome
        @test new_Ψ.args[index].outcome_extra_covariates == (Symbol(22001),)
        @test new_Ψ.args[index].treatment_confounders == (RSID_198 = (:PC1, :PC2), RSID_2 = (:PC1, :PC2),)
    end
    @test new_Ψ.args[1].treatment_values == (RSID_198 = "GA", RSID_2 = "GG")
    @test new_Ψ.args[2].treatment_values == (RSID_198 = "GA", RSID_2 = "AA")
    
    # If the allele is not present 
    variants_alleles = Dict(:RSID_198 => Set(["AA"]))
    @test_throws TargeneCore.AbsentAlleleError("RSID_198", "AG") TargeneCore.adjust_estimand_fields(Ψ, variants_alleles, pcs)
end

#####################################################################
###############           END-TO-END TESTS            ###############
#####################################################################

@testset "Test inputs_from_estimands" begin
    # Genotypes encoded as strings
    # No batching of estimands files
    # No positivity constraint
    tmpdir = mktempdir()
    estimands_filename = make_estimands_configuration_file()
    copy!(ARGS, [
        "estimation-inputs",
        estimands_filename,
        string("--traits-file=", joinpath(TESTDIR, "data", "traits_1.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--call-threshold=0.8",  
        "--verbosity=0",
        "--positivity-constraint=0"
    ])
    TargeneCore.julia_main()
    ## Data File
    data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test names(data) == [
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", 
        "CONTINUOUS_2", "COV_1", "21003", "22001", "TREAT_1", 
        "PC1", "PC2", "RSID_2", "RSID_198"
    ]
    @test size(data) == (490, 13)

    ## Estimands file:
    output_estimands = deserialize(joinpath(tmpdir, "final.estimands.jls")).estimands
    # There are 5 initial estimands containing a :ALL
    # Those are duplicated for each of the 4 outcomes.
    @test length(output_estimands) == 20
    # In all cases the PCs are appended to the confounders.
    for Ψ ∈ output_estimands
        # Input Estimand 1
        if Ψ isa TMLE.StatisticalIATE
            @test Ψ.treatment_confounders == (RSID_2 = (:PC1, :PC2), TREAT_1 = (:PC1, :PC2))
        
        # Input Estimand 2
        elseif Ψ isa TMLE.StatisticalATE && Ψ.treatment_values == (RSID_2 = (case = "AA", control = "GG"),)
            @test Ψ.treatment_confounders == (RSID_2 = (Symbol("22001"), :PC1, :PC2),)
            @test Ψ.outcome_extra_covariates == (Symbol("21003"), Symbol("COV_1"))

        # Input Estimand 3
        elseif Ψ isa TMLE.StatisticalCM && Ψ.treatment_values == (RSID_2 = "AA", )
            @test Ψ.treatment_confounders == (RSID_2 = (Symbol("22001"), :PC1, :PC2),)
            @test Ψ.outcome_extra_covariates == (Symbol("21003"), Symbol("COV_1"))

        # Input Estimand 4
        elseif Ψ isa TMLE.StatisticalATE && Ψ.treatment_values == (RSID_198 = (case = "AG", control = "AA"), RSID_2 = (case = "AA", control = "GG"))
            @test Ψ.treatment_confounders == (RSID_198 = (:PC1, :PC2), RSID_2 = (:PC1, :PC2))
            @test Ψ.outcome_extra_covariates == (Symbol("22001"),)
        
        # Input Estimand 5: GA is corrected to AG to match the data
        elseif Ψ isa JointEstimand
            @test Ψ.args[1].treatment_values == (RSID_198 = "AG", RSID_2 = "GG")
            @test Ψ.args[2].treatment_values == (RSID_198 = "AG", RSID_2 = "AA")
            @test Ψ.args[1].treatment_confounders == Ψ.args[2].treatment_confounders == (RSID_198 = (:PC1, :PC2), RSID_2 = (:PC1, :PC2))
            @test Ψ.args[1].outcome_extra_covariates == Ψ.args[2].outcome_extra_covariates == (Symbol("22001"),)
        else
            throw(AssertionError(string("Which input did this output come from: ", Ψ)))
        end
    end

    η_counts = TMLE.nuisance_function_counts(output_estimands)
    memcost, _ = TMLE.evaluate_proxy_costs(output_estimands, η_counts)
    shuffled_estimands = shuffle(MersenneTwister(123), output_estimands)
    shuffled_memcost, _ = TMLE.evaluate_proxy_costs(shuffled_estimands, η_counts)
    @test memcost < shuffled_memcost

    # Increase positivity constraint
    copy!(ARGS, [
        "estimation-inputs",
        estimands_filename,
        string("--traits-file=", joinpath(TESTDIR, "data", "traits_1.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--call-threshold=0.8",  
        "--verbosity=0",
        "--positivity-constraint=0.01"
    ])
    TargeneCore.julia_main()

    # The IATES are the most sensitives
    outestimands = deserialize(joinpath(tmpdir, "final.estimands.jls")).estimands
    @test all(Ψ isa Union{TMLE.StatisticalCM, TMLE.StatisticalATE, JointEstimand} for Ψ in outestimands)
    @test size(outestimands, 1) == 16

    # No Remaining Estimand Error
    copy!(ARGS, [
        "estimation-inputs",
        estimands_filename,
        string("--traits-file=", joinpath(TESTDIR, "data", "traits_1.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--call-threshold=0.8",  
        "--verbosity=0",
        "--positivity-constraint=1"
    ])
    
    @test_throws TargeneCore.NoRemainingParamsError(1.) TargeneCore.julia_main()
end

@testset "Test estimation_inputs from-estimands-file: no wildcard" begin
    estimands_filename = make_estimands_configuration_file(make_estimands_configuration_no_wildcard)
    tmpdir = mktempdir()
    copy!(ARGS, [
        "estimation-inputs",
        estimands_filename,
        string("--traits-file=", joinpath(TESTDIR, "data", "traits_1.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb")),
        string("--outprefix=", joinpath(tmpdir, "final")),
        "--batchsize=2",
        "--verbosity=0",
        "--positivity-constraint=0."
    ])
    TargeneCore.julia_main()

    ## Data File
    data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test names(data) == [
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", 
        "CONTINUOUS_2", "COV_1", "21003", "22001", "TREAT_1", 
        "PC1", "PC2", "RSID_2"
    ]
    @test size(data) == (490, 12)

    ## Parameter files:
    # There are 3 initial estimands, 1 containing a *
    # that will be duplicated for each of the 4 targets.
    # The PCs are appended to the confounders.
    # Estimands are batched by 2
    output_estimands = [
        deserialize(joinpath(tmpdir, "final.estimands_1.jls")).estimands,
        deserialize(joinpath(tmpdir,"final.estimands_1.jls")).estimands,
        deserialize(joinpath(tmpdir,"final.estimands_1.jls")).estimands
    ]
    for index in 1:3
        @test length(output_estimands[index]) == 2
    end
end


end

true