module TestLocoGwasEstimands

using Test
using SnpArrays
using TargeneCore
using Arrow
using DataFrames
using Serialization
using TMLE
using CSV

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "tl_inputs", "test_utils.jl"))

function get_summary_stats(estimands)
    outcomes = [only(TargeneCore.outcome_variables(Ψ)) for Ψ in estimands]
    results = DataFrame(ESTIMAND = estimands, OUTCOME = outcomes)
    return sort(combine(groupby(results, :OUTCOME), nrow), :OUTCOME)
end

function check_estimands_levels_order(estimands)
    for Ψ in estimands
        # If the two components are present, the first is the 0 -> 1 and the second is the 1 -> 2
        if length(Ψ.args) == 2
            @test Ψ.args[1].treatment_values[1] == (case = 0x01, control = 0x00)
            @test Ψ.args[2].treatment_values[1] == (case = 0x02, control = 0x01)
        else
            # Otherwise we check they are one or the other
            arg = only(Ψ.args)
            @test arg.treatment_values[1]==(case = 0x01, control = 0x00) ||
            arg.treatment_values[1]==(case = 0x02, control = 0x01)
        end
   end
end
@testset "Test loco-gwas serial: no positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 5,
        "positivity-constraint" => 0.0,

        "%COMMAND%" => "allele-independent", 

        "allele-independent" => Dict{String, Any}(
            "call-threshold" => nothing,
            "config" => joinpath(TESTDIR, "data", "config_gwas_serial.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "ukbb_traits.csv"),
            "pcs" => joinpath(TESTDIR, "data", "ukbb_pcs.csv"),
            "genotype-prefix" => joinpath(TESTDIR, "data", "ukbb", "genotypes" ,"ukbb_1."), 
        ),
    )
    tl_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(trait_data) == (1940, 886)
    
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa JointEstimand for e in estimands)

    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:BINARY_1, :BINARY_2, :CONTINUOUS_1, :CONTINUOUS_2, :TREAT_1], 
        nrow = repeat([875], 5)
    )

    check_estimands_levels_order(estimands)
end

@testset "Test loco-gwas serial: positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 5,
        "positivity-constraint" => 0.2,

        "%COMMAND%" => "allele-independent", 

        "allele-independent" => Dict{String, Any}(
            "call-threshold" => nothing,
            "config" => joinpath(TESTDIR, "data", "config_gwas_serial.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "ukbb_traits.csv"),
            "pcs" => joinpath(TESTDIR, "data", "ukbb_pcs.csv"),
            "genotype-prefix" => joinpath(TESTDIR, "data", "ukbb", "genotypes" ,"ukbb_1")
        ),
    )
    tl_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(trait_data) == (1940, 886)
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa JointEstimand for e in estimands)
    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:BINARY_1, :BINARY_2, :CONTINUOUS_1, :CONTINUOUS_2, :TREAT_1], 
        nrow = repeat([777], 5)
    )
   
    check_estimands_levels_order(estimands)

end

@testset "Test loco-gwas parallel: no positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 5,
        "positivity-constraint" => 0.0,

        "%COMMAND%" => "allele-independent", 

        "allele-independent" => Dict{String, Any}(
            "call-threshold" => nothing,
            "config" => joinpath(TESTDIR, "data", "config_gwas_parallel.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "ukbb_traits.csv"),
            "pcs" => joinpath(TESTDIR, "data", "ukbb_pcs.csv"),
            "genotype-prefix" => joinpath(TESTDIR, "data", "ukbb", "genotypes" ,"ukbb_1."), 
        ),
    )
    tl_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(trait_data) == (1940, 886)
    
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa JointEstimand for e in estimands)

    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:BINARY_1, :BINARY_2, :CONTINUOUS_1, :CONTINUOUS_2, :TREAT_1], 
        nrow = repeat([875], 5)
    )

    check_estimands_levels_order(estimands)
end

@testset "Test loco-gwas parallel: positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 5,
        "positivity-constraint" => 0.2,

        "%COMMAND%" => "allele-independent", 

        "allele-independent" => Dict{String, Any}(
            "call-threshold" => nothing,
            "config" => joinpath(TESTDIR, "data", "config_gwas_parallel.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "ukbb_traits.csv"),
            "pcs" => joinpath(TESTDIR, "data", "ukbb_pcs.csv"),
            "genotype-prefix" => joinpath(TESTDIR, "data", "ukbb", "genotypes" ,"ukbb_1")
        ),
    )
    tl_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(trait_data) == (1940, 886)
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa JointEstimand for e in estimands)
    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:BINARY_1, :BINARY_2, :CONTINUOUS_1, :CONTINUOUS_2, :TREAT_1], 
        nrow = repeat([777], 5)
    )
   
    check_estimands_levels_order(estimands)

end

end
true