module TestPermutationTest

using Test
using CSV
using DataFrames
using StableRNGs
using TargeneCore
using TMLE
using Arrow
using Serialization

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "testutils.jl"))

@testset "Test permuted_estimand!" begin
    estimates = make_estimates()
    Ψ = estimates[1].TMLE.estimand
    @test Ψ isa TMLE.StatisticalIATE
    # Treatment and Outcome
    permutation_variables = Set([Ψ.outcome, :RSID_103])
    Ψpermuted = TargeneCore.permuted_estimand!(permutation_variables, Ψ)
    @test Ψpermuted.outcome == TargeneCore.permuted_name(Ψ.outcome)
    @test keys(Ψpermuted.treatment_values) == (:RSID_103_permuted, :rs10043934)
    @test values(Ψpermuted.treatment_values) == values(Ψ.treatment_values)
    @test keys(Ψpermuted.treatment_confounders) == (:RSID_103_permuted, :rs10043934)
    @test values(Ψpermuted.treatment_confounders) == values(Ψ.treatment_confounders)
    @test Ψpermuted.outcome_extra_covariates == Ψ.outcome_extra_covariates
    # Treatment only
    permutation_variables = Set([:rs10043934])
    Ψpermuted = TargeneCore.permuted_estimand!(permutation_variables, Ψ)
    @test Ψpermuted.outcome == Ψ.outcome
    @test keys(Ψpermuted.treatment_values) == (:RSID_103, :rs10043934_permuted)
    @test values(Ψpermuted.treatment_values) == values(Ψ.treatment_values)
    @test keys(Ψpermuted.treatment_confounders) == (:RSID_103, :rs10043934_permuted)
    @test values(Ψpermuted.treatment_confounders) == values(Ψ.treatment_confounders)
    @test Ψpermuted.outcome_extra_covariates == Ψ.outcome_extra_covariates
    # Outcome only
    permutation_variables = Set([Ψ.outcome])
    Ψpermuted = TargeneCore.permuted_estimand!(permutation_variables, Ψ)
    @test Ψpermuted.outcome == TargeneCore.permuted_name(Ψ.outcome)
    @test Ψpermuted.treatment_values == Ψ.treatment_values
    @test Ψpermuted.treatment_confounders == Ψ.treatment_confounders
    @test Ψpermuted.outcome_extra_covariates == Ψ.outcome_extra_covariates
    # Composed Estimand
    Ψ = estimates[3].TMLE.estimand
    @test Ψ isa JointEstimand
    outcome = Symbol("High light scatter reticulocyte percentage")
    permutation_variables = Set([:RSID_103, outcome])
    Ψpermuted = TargeneCore.permuted_estimand!(permutation_variables, Ψ)
    arg₁ = Ψpermuted.args[1]
    arg₂ = Ψpermuted.args[2]
    @test arg₁.outcome == arg₂.outcome == TargeneCore.permuted_name(outcome)
    @test keys(arg₁.treatment_values) == keys(arg₂.treatment_values) == (:RSID_103_permuted, :rs10043934)
    @test values(arg₁.treatment_values) == values(Ψ.args[1].treatment_values)
    @test values(arg₂.treatment_values) == values(Ψ.args[2].treatment_values)
    @test arg₂.outcome_extra_covariates == arg₁.outcome_extra_covariates == Ψ.args[1].outcome_extra_covariates
end

@testset "Test make_permutation_parameters" begin
    estimands = [Ψ.TMLE.estimand for Ψ ∈ make_estimates()[1:end-1]]
    expected_permuted_variables = Set([
        :rs117913124,
        :rs10043934,
        :RSID_104,
        Symbol("L50-L54 Urticaria and erythema"),
        Symbol("High light scatter reticulocyte percentage"),
        :RSID_103
    ])
    # Each estimand has 2 treatment variables and 1 outcome
    # Each estimand will thus lead to 3 new permutation estimands
    # A total of 12
    permutation_estimands, all_permuted_variables = TargeneCore.make_permutation_parameters(estimands; optimize=false, orders=(1,))
    @test length(permutation_estimands) == 12
    @test all(values(Ψ.treatment_values) == values(estimands[1].treatment_values) for Ψ ∈ permutation_estimands[1:3])
    @test all(values(Ψ.treatment_values) == values(estimands[2].treatment_values) for Ψ ∈ permutation_estimands[4:6])
    @test all(Ψ isa JointEstimand for Ψ ∈ permutation_estimands[7:9])
    @test all(values(Ψ.treatment_values) == values(estimands[4].treatment_values) for Ψ ∈ permutation_estimands[10:12])

    @test all_permuted_variables == expected_permuted_variables
    # For each of the 4 parameters there are:
    # - 3 x order 1 permutations
    # - 3 x order 2 permutations
    # - 1 order 3 permutation
    # = 7 parameters
    # Total: 7*4 = 28 permutation estimands
    permutation_estimands, all_permuted_variables = TargeneCore.make_permutation_parameters(estimands; optimize=false, orders=(1, 2, 3))
    @test length(permutation_estimands) == 28
    @test all(values(Ψ.treatment_values) == values(estimands[1].treatment_values) for Ψ ∈ permutation_estimands[1:7])
    @test all(values(Ψ.treatment_values) == values(estimands[2].treatment_values) for Ψ ∈ permutation_estimands[8:14])
    @test all(Ψ isa JointEstimand for Ψ ∈ permutation_estimands[15:21])
    @test all(values(Ψ.treatment_values) == values(estimands[4].treatment_values) for Ψ ∈ permutation_estimands[22:28])

    @test all_permuted_variables == expected_permuted_variables
end

@testset "Test filter_by_positivity_threshold" begin
    dataset = make_dataset()
    estimands = make_estimands()
    permuted_estimands, permuted_variables = TargeneCore.make_permutation_parameters(estimands)
    # No positivity constraint
    valid_estimands = TargeneCore.permute_dataset_and_get_valid_estimands!(
        dataset, 
        permuted_variables, 
        permuted_estimands; 
        rng_seed=123, 
        max_attempts=1,
        positivity_constraint=nothing
    )
    @test Set(names(dataset)) == Set([
        "rs10043934",
        "rs117913124",
        "RSID_103",
        "RSID_104",
        "High light scatter reticulocyte percentage",
        "L50-L54 Urticaria and erythema",
        "rs117913124_permuted",
        "rs10043934_permuted",
        "RSID_104_permuted",
        "L50-L54 Urticaria and erythema_permuted",
        "High light scatter reticulocyte percentage_permuted",
        "RSID_103_permuted"
    ])
    @test valid_estimands == permuted_estimands
    # Positivity constraint, 1 attempt
    dataset = make_dataset()
    valid_estimands = TargeneCore.permute_dataset_and_get_valid_estimands!(
        dataset, 
        permuted_variables, 
        permuted_estimands; 
        rng_seed=123, 
        max_attempts=1,
        positivity_constraint=0.05
    )
    @test length(valid_estimands) < length(permuted_estimands)
    # Positivity constraint, 5 attempt
    dataset = make_dataset()
    rng = StableRNG(123)
    log_sequence = [
        (:info, "Initial permutation resulted in a loss of estimands, attempting again."),
        (:info, "After attempt 2/5 #Estimands=8/12"),
        (:info, "After attempt 3/5 #Estimands=8/12"),
        (:info, "After attempt 4/5 #Estimands=8/12"),
        (:info, "After attempt 5/5 #Estimands=9/12"),
        ]
    valid_estimands = @test_logs log_sequence... TargeneCore.permute_dataset_and_get_valid_estimands!(
        dataset, 
        permuted_variables, 
        permuted_estimands; 
        rng_seed=123, 
        max_attempts=5,
        positivity_constraint=0.05,
        verbosity=1
    );
    @test length(valid_estimands) == 9
end

@testset "Test run_permutation_test: no positivity constraint" begin
    make_fake_outputs()
    parsed_args = Dict(
        "batch-size" => 5,
        "verbosity" => 0,
        "positivity-constraint" => nothing,
        "out-prefix" => ".",

        "%COMMAND%" => "permutation-tests",

        "permutation-tests" => Dict{String, Any}(
            "dataset" => joinpath(TESTDIR, "data", "dataset.csv"),
            "results" => "tmle_output.hdf5",
            "estimator-key" => "TMLE",
            "pval-threshold" => 1e-10,
            "limit" => nothing,
            "rng" => 123,
            "orders" => "1,2",
            "max-attempts" => 1
        ),
    )
    tl_inputs(parsed_args)

    # Check permutation dataset file
    dataset_path = TargeneCore.filepath_from_prefix("."; filename="permutation_dataset.arrow")
    data = DataFrame(Arrow.Table(dataset_path))
    permuted_cols = filter(x -> endswith(x, "permuted"), names(data))
    @test Set(permuted_cols) == Set([
        "rs10043934_permuted",
        "High light scatter reticulocyte percentage_permuted",
        "RSID_103_permuted"
    ])
    
    # Only two estimates pass the threshold
    # Each estimate produces 6 new estimands
    # Total: 12 split into 3 files of chunksize 5
    batch_1 = deserialize("permutation_estimands_1.jls").estimands
    @test length(batch_1) == 5
    batch_2 = deserialize("permutation_estimands_2.jls").estimands
    @test length(batch_2) == 5
    batch_3 = deserialize("permutation_estimands_3.jls").estimands
    @test length(batch_3) == 2

    # Clean
    clean()
    rm("permutation_estimands_1.jls")
    rm("permutation_estimands_2.jls")
    rm("permutation_estimands_3.jls")
    rm("permutation_dataset.arrow")
end

@testset "Test run_permutation_test: positivity constraint" begin
    make_fake_outputs()
    parsed_args = Dict(
        "batch-size" => 50,
        "verbosity" => 0,
        "positivity-constraint" => 0.5,
        "out-prefix" => ".",

        "%COMMAND%" => "permutation-tests",

        "permutation-tests" => Dict{String, Any}(
            "dataset" => joinpath(TESTDIR, "data", "dataset.csv"),
            "results" => "tmle_output.hdf5",
            "estimator-key" => "TMLE",
            "pval-threshold" => 1e-10,
            "limit" => nothing,
            "rng" => 123,
            "orders" => "1,2",
            "max-attempts" => 5
        ),
    )
    # None passs the positivity constraint => throws error
    @test_throws ErrorException("No permuted estimand remaining, consider increasing the p-value threshold or the maximum number of attempts.") tl_inputs(parsed_args)

    parsed_args["positivity-constraint"] = 0.1
    tl_inputs(parsed_args)
    # Check permutation dataset file
    dataset_path = TargeneCore.filepath_from_prefix("."; filename="permutation_dataset.arrow")
    data = DataFrame(Arrow.Table(dataset_path))
    permuted_cols = filter(x -> endswith(x, "permuted"), names(data))
    @test Set(permuted_cols) == Set([
        "rs10043934_permuted",
        "High light scatter reticulocyte percentage_permuted",
        "RSID_103_permuted"
    ])
    
    # Only two estimates pass the pvalue threshold, so a maximum of 12
    # and then only 4 make it due to positivity constraint
    estimands = deserialize("permutation_estimands_1.jls").estimands
    @test length(estimands) < 12

    # Clean
    clean()
    rm("permutation_estimands_1.jls")
    rm("permutation_dataset.arrow")
end

end

true