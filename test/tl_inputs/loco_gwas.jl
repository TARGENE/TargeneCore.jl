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

function summary_stats_df(estimands)
    estimands = DataFrame(ESTIMAND=estimands)
    estimands.ESTIMAND_TYPE = [string(typeof(first(x.args))) for x in estimands.ESTIMAND]
    estimands.ORDER = [length(first(x.args).treatment_values) for x in estimands.ESTIMAND]
    return sort(DataFrames.combine(groupby(estimands, [:ESTIMAND_TYPE, :ORDER]), nrow), [:ORDER, :ESTIMAND_TYPE])
end
@testset "Test loco-gwas from flat list: no positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 5,
        "positivity-constraint" => 0.0,

        "%COMMAND%" => "loco-gwas", 

        "loco-gwas" => Dict{String, Any}(
            "config" => joinpath(TESTDIR, "data", "config_gwas.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "ukbb_traits.csv"),
            "pcs" => joinpath(TESTDIR, "data", "ukbb_pcs.csv"),
            "bed-file" => joinpath(TESTDIR, "data", "ukbb", "genotypes" ,"ukbb_1.bed"), 
            "bim-file" => joinpath(TESTDIR, "data", "ukbb", "genotypes" ,"ukbb_1.bim"),
            "fam-file" => joinpath(TESTDIR, "data", "ukbb", "genotypes" ,"ukbb_1.fam")
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
    summary_stats = summary_stats_df(estimands)
    @test summary_stats == DataFrame(
        ESTIMAND_TYPE=["TMLE.StatisticalATE"],
        ORDER=[1],
        nrow=[868]
    )
   @test all((case = 0x01, control = 0x00) == e.args[1].treatment_values[1] for e in estimands)
   @test all((case = 0x02, control = 0x01) == e.args[2].treatment_values[1] for e in estimands)
        
    
end

@testset "Test loco-gwas from flat list: positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 5,
        "positivity-constraint" => 0.2,

        "%COMMAND%" => "loco-gwas", 

        "loco-gwas" => Dict{String, Any}(
            "config" => joinpath(TESTDIR, "data", "config_gwas.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "ukbb_traits.csv"),
            "pcs" => joinpath(TESTDIR, "data", "ukbb_pcs.csv"),
            "bed-file" => joinpath(TESTDIR, "data", "ukbb", "genotypes" ,"ukbb_1.bed"), 
            "bim-file" => joinpath(TESTDIR, "data", "ukbb", "genotypes" ,"ukbb_1.bim"),
            "fam-file" => joinpath(TESTDIR, "data", "ukbb", "genotypes" ,"ukbb_1.fam")
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
    summary_stats = summary_stats_df(estimands)
    @test summary_stats == DataFrame(
        ESTIMAND_TYPE=["TMLE.StatisticalATE"],
        ORDER=[1],
        nrow=[777]
    )
   
    # because some estimands can be dropped here check if all estimands 
    # match the proper form at the first index
    @test all(e.args[1].treatment_values[1]==(case = 0x01, control = 0x00) ||
    e.args[1].treatment_values[1]==(case = 0x02, control = 0x01) for e in estimands)
    
    # if both estimands survive positivity check if the second estimand is
    # of the form Aa -> AA
    survived_estimands = []

    for e in estimands
        if length(e.args) == 2
            push!(survived_estimands, e)
        end
    end
    @test all((case = 0x02, control = 0x01) == e.args[2].treatment_values[1] for e in survived_estimands)

end

end
true