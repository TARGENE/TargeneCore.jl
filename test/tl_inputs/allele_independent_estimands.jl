module TestAlleleIndependeEstimands

using Test
using TargeneCore
using Arrow
using DataFrames
using Serialization
using TMLE

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "tl_inputs", "test_utils.jl"))

function summary_stats_df(estimands)
    estimands = DataFrame(ESTIMAND=estimands)
    estimands.ESTIMAND_TYPE = [string(typeof(first(x.args))) for x in estimands.ESTIMAND]
    estimands.ORDER = [length(first(x.args).treatment_values) for x in estimands.ESTIMAND]
    return sort(DataFrames.combine(groupby(estimands, [:ESTIMAND_TYPE, :ORDER]), nrow), [:ORDER, :ESTIMAND_TYPE])
end
@testset "Test treatment_tuples_from_groups" begin
    treatments_list = [
        [:RSID_1, :RSID_2],
        [:RSID_3, :RSID_4],
        [:RSID_5],
        ]
    order_1 = TargeneCore.treatment_tuples_from_groups(treatments_list, [1])
    @test order_1 == [
        (:RSID_1,),
        (:RSID_2,),
        (:RSID_3,),
        (:RSID_4,),
        (:RSID_5,)
    ]
    order_2 = TargeneCore.treatment_tuples_from_groups(treatments_list, [2])
    @test order_2 == [
        (:RSID_1, :RSID_3),
        (:RSID_1, :RSID_4),
        (:RSID_1, :RSID_5),
        (:RSID_2, :RSID_3),
        (:RSID_2, :RSID_4),
        (:RSID_2, :RSID_5),
        (:RSID_3, :RSID_5),
        (:RSID_4, :RSID_5)
    ]
    order_3 = TargeneCore.treatment_tuples_from_groups(treatments_list, [3])
    @test order_3 == [
        (:RSID_1, :RSID_3, :RSID_5),
        (:RSID_1, :RSID_4, :RSID_5),
        (:RSID_2, :RSID_3, :RSID_5),
        (:RSID_2, :RSID_4, :RSID_5)
    ]
    order_2_3 = TargeneCore.treatment_tuples_from_groups(treatments_list, [2, 3])
    @test order_2_3 == sort(vcat(order_2, order_3))
end

@testset "Test misc" begin
    @test TargeneCore.default_order(ATE) == TargeneCore.default_order(CM) == [1]
    @test TargeneCore.default_order(IATE) == [2]

    treatment_tuples = TargeneCore.treatment_tuples_from_single_list([:T_1, :T_2], [1])
    @test treatment_tuples == [[:T_1], [:T_2]]
    treatment_tuples = TargeneCore.treatment_tuples_from_single_list([:T_1, :T_2, :T_3], [2, 3])
    @test treatment_tuples == [
        [:T_1, :T_2],
        [:T_1, :T_3],
        [:T_2, :T_3],
        [:T_1, :T_2, :T_3]
    ]
end

@testset "Test allele-independent from groups: no positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 100,
        "verbosity" => 0,
        "positivity-constraint" => nothing,

        "%COMMAND%" => "allele-independent", 

        "allele-independent" => Dict{String, Any}(
            "config" => joinpath(TESTDIR, "data", "config_groups.yaml"),
            "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
            "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
            "call-threshold" => 0.8,  
            "bgen-prefix" => joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb"), 
            ),
    )
    tl_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test sort(names(trait_data)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2", "RSID_2", "RSID_102", 
        "RSID_17", "RSID_198", "RSID_99"])
    @test size(trait_data) == (490, 16)
    # Check estimands
    estimands = deserialize(joinpath(tmpdir, "final.estimands_1.jls")).estimands
    summary_stats = summary_stats_df(estimands)
    @test summary_stats == DataFrame(
        ESTIMAND_TYPE = ["TMLE.StatisticalCM", "TMLE.StatisticalIATE", "TMLE.StatisticalIATE"],
        ORDER = [2, 2, 3],
        nrow = [40, 40, 16]
    )
    # The total number of estimands is 56
    n_traits = 4
    n_treat_comb_order_2 = 5
    n_treat_comb_order_3 = 2
    ngroups = 2
    @test ngroups*n_traits*(2n_treat_comb_order_2+n_treat_comb_order_3) == 96
end

@testset "Test allele-independent from groups: with positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 10,
        "positivity-constraint" => 0.01,

        "%COMMAND%" => "allele-independent", 

        "allele-independent" => Dict{String, Any}(
            "config" => joinpath(TESTDIR, "data", "config_groups.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
            "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
            "call-threshold" => 0.8,  
            "bgen-prefix" => joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb"), 
        ),
    )
    tl_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test sort(names(trait_data)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2", "RSID_2", "RSID_102", 
        "RSID_17", "RSID_198", "RSID_99"])
    @test size(trait_data) == (490, 16)
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa ComposedEstimand for e in estimands)

    summary_stats = summary_stats_df(estimands)
    @test summary_stats == DataFrame(
        ESTIMAND_TYPE=["TMLE.StatisticalCM", "TMLE.StatisticalIATE", "TMLE.StatisticalIATE"],
        ORDER = [2, 2, 3],
        nrow=[40, 40, 8]
    )
end

@testset "Test allele-independent from flat list: no positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 5,
        "positivity-constraint" => 0.0,

        "%COMMAND%" => "allele-independent", 

        "allele-independent" => Dict{String, Any}(
            "config" => joinpath(TESTDIR, "data", "config_flat.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
            "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
            "call-threshold" => 0.8,  
            "bgen-prefix" => joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb"), 
        ),
    )
    tl_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test sort(names(trait_data)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2",
        "RSID_102", "RSID_99", "RSID_17"])
    @test size(trait_data) == (490, 14)
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa ComposedEstimand for e in estimands)
    summary_stats = summary_stats_df(estimands)
    @test summary_stats == DataFrame(
        ESTIMAND_TYPE=["TMLE.StatisticalATE", "TMLE.StatisticalCM", "TMLE.StatisticalIATE", "TMLE.StatisticalIATE"],
        ORDER=[1, 1, 2, 3],
        nrow=[16, 16, 24, 16]
    )
end

@testset "Test allele-independent from flat list: positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => nothing,
        "positivity-constraint" => 0.01,

        "%COMMAND%" => "allele-independent", 

        "allele-independent" => Dict{String, Any}(
            "config" => joinpath(TESTDIR, "data", "config_flat.yaml"), 
            "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
            "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
            "call-threshold" => 0.8,  
            "bgen-prefix" => joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb"), 
        ),
    )
    tl_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test sort(names(trait_data)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2",
        "RSID_102", "RSID_99", "RSID_17"])
    @test size(trait_data) == (490, 14)
    # Check estimands
    estimands = deserialize(joinpath(tmpdir, "final.estimands_1.jls")).estimands
    @test all(e isa ComposedEstimand for e in estimands)
    summary_stats = summary_stats_df(estimands)
    @test summary_stats == DataFrame(
        ESTIMAND_TYPE=["TMLE.StatisticalATE", "TMLE.StatisticalCM", "TMLE.StatisticalIATE", "TMLE.StatisticalIATE"],
        ORDER=[1, 1, 2, 3],
        nrow=[16, 16, 24, 4]
        )
end

end

true