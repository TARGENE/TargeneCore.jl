module TestGwisEstimands

using Test
using SnpArrays
using TargeneCore
using Arrow
using DataFrames
using Serialization
using TMLE
using CSV
using YAML
using StatsBase

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "testutils.jl"))

function get_summary_stats(estimands)
    outcomes = [TargeneCore.get_outcome(Ψ) for Ψ in estimands]
    results = DataFrame(ESTIMAND = estimands, OUTCOME = outcomes)
    return sort(combine(groupby(results, :OUTCOME), nrow), :OUTCOME)
end

function check_estimands_levels_interactions(estimands)
    string_treatments = YAML.load_file(joinpath(TESTDIR, "data", "config_gwis.yaml"))["extra_treatments"]
    extra_treatments = Set(Symbol.(string_treatments))
    ATE_count, GxE_Count, GxG_Count, higher_order_count = 0, 0, 0, 0

    for Ψ in estimands
        treatment_set = Set(keys(Ψ.args[1].treatment_values))
        # Find the variant (treatment that's not in extra_treatments)
        variants = setdiff(treatment_set, extra_treatments)
        # Will crash if there are multiple variants in the set
        # If the variant is used for interaction e.g length(variants) == 0 after the setdiff
        # Account for this by just picking the first treatment in tests to ensure that single effect of said variant works correctly
        if isempty(variants)
            variant = first(treatment_set)
        else
            variant = only(variants)
        end
        
        # For all estimands, just check that each component has a valid variant transition
        for arg in Ψ.args
            @test arg.treatment_values[variant] == (control = 0x00, case = 0x01) ||
                  arg.treatment_values[variant] == (control = 0x01, case = 0x02)
        end
        
        if Ψ.args[1] isa TMLE.StatisticalATE
            @test length(Ψ.args[1].treatment_values) == 1
            ATE_count += 1
        elseif haskey(Ψ.args[1].treatment_values, :TREAT_1) && length(Ψ.args[1].treatment_values) == 2
            GxE_Count += 1
        elseif haskey(Ψ.args[1].treatment_values, :rs3693265) && length(Ψ.args[1].treatment_values) == 2
            GxG_Count += 1
        elseif haskey(Ψ.args[1].treatment_values, :rs3693265) && haskey(Ψ.args[1].treatment_values, :TREAT_1) && length(Ψ.args[1].treatment_values) == 3
            higher_order_count += 1
        else
            @test false
        end
    end
end

@testset "Test inputs_from_config gwis: no positivity constraint" begin
    tmpdir = mktempdir()
    copy!(ARGS, [
        "estimation-inputs",
        joinpath(TESTDIR, "data", "config_gwis.yaml"),
        string("--traits-file=", joinpath(TESTDIR, "data", "ukbb_traits.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "ukbb_pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "genotypes" , "ukbb_1.")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--batchsize=5",
        "--verbosity=1",
        "--positivity-constraint=0"
    ])
    TargeneCore.julia_main()
    # Check dataset
    dataset = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(dataset) == (1940, 886)
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa JointEstimand for e in estimands)

    summary_stats = get_summary_stats(estimands)
    # Removes estimands that have no representation in the data as well as redundant interaction (e.g. variant treatment and itself)
    @test summary_stats == DataFrame(
        OUTCOME = [:BINARY_1, :BINARY_2, :CONTINUOUS_1, :CONTINUOUS_2], 
        nrow = repeat([3493], 4)
    )
    check_estimands_levels_interactions(estimands)
end


@testset "Test inputs_from_config gwis: positivity constraint" begin
    tmpdir = mktempdir()
    copy!(ARGS, [
        "estimation-inputs",
        joinpath(TESTDIR, "data", "config_gwis.yaml"),
        string("--traits-file=", joinpath(TESTDIR, "data", "ukbb_traits.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "ukbb_pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "genotypes" , "ukbb_1.")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--batchsize=5",
        "--verbosity=0",
        "--positivity-constraint=0.2"
    ])
    TargeneCore.julia_main()
    # Check dataset
    dataset = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(dataset) == (1940, 886)
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    # The positivity constraint reduces the number of variants
    @test all(e isa JointEstimand for e in estimands)
    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:BINARY_1, :BINARY_2, :CONTINUOUS_1, :CONTINUOUS_2], 
        nrow = repeat([923], 4)
    )
   
    check_estimands_levels_interactions(estimands)
end

@testset "Test inputs_from_config gwis: incorrect estimand type" begin
    tmpdir = mktempdir()
    copy!(ARGS, [
        "estimation-inputs",
        joinpath(TESTDIR, "data", "config_gwis_problematic.yaml"),
        string("--traits-file=", joinpath(TESTDIR, "data", "ukbb_traits.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "ukbb_pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "genotypes" , "ukbb_1.")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--batchsize=5",
        "--verbosity=0",
        "--positivity-constraint=0"
    ])
    @test_throws TaskFailedException TargeneCore.julia_main()
end

end

true