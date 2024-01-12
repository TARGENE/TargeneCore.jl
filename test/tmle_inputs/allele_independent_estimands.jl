module TestAlleleIndependeEstimands

using Test
using TargeneCore
using Arrow
using DataFrames
using Serialization

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "tmle_inputs", "test_utils.jl"))

@testset "Test generate_treatments_combinations" begin
    treatments_list = [
        [:RSID_1, :RSID_2],
        [:RSID_3, :RSID_4],
        [:RSID_5],
        ]
    order_1 = TargeneCore.generate_treatments_combinations(treatments_list, [1])
    @test order_1 == [
        (:RSID_1,),
        (:RSID_2,),
        (:RSID_3,),
        (:RSID_4,),
        (:RSID_5,)
    ]
    order_2 = TargeneCore.generate_treatments_combinations(treatments_list, [2])
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
    order_3 = TargeneCore.generate_treatments_combinations(treatments_list, [3])
    @test order_3 == [
        (:RSID_1, :RSID_3, :RSID_5),
        (:RSID_1, :RSID_4, :RSID_5),
        (:RSID_2, :RSID_3, :RSID_5),
        (:RSID_2, :RSID_4, :RSID_5)
    ]
    order_2_3 = TargeneCore.generate_treatments_combinations(treatments_list, [2, 3])
    @test order_2_3 == sort(vcat(order_2, order_3))
end

@testset "Test allele-independent: no positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "allele-independent" => Dict{String, Any}(
            "config" => joinpath(TESTDIR, "data", "interaction_config.yaml"), 
            ),
        "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
        "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
        "call-threshold" => 0.8,  
        "%COMMAND%" => "allele-independent", 
        "bgen-prefix" => joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => 100,
        "positivity-constraint" => nothing,
    )
    tmle_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test sort(names(trait_data)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2", "RSID_2", "RSID_102", 
        "RSID_17", "RSID_198", "RSID_99"])
    @test size(trait_data) == (490, 16)
    # Check estimands
    tf_estimands = Dict(:TF1 => [], :TF2 => [])
    for file in readdir(tmpdir)
        if startswith(file, "final.TF1")
            append!(tf_estimands[:TF1], deserialize(joinpath(tmpdir, file)).estimands)
        elseif startswith(file, "final.TF2")
            append!(tf_estimands[:TF2], deserialize(joinpath(tmpdir, file)).estimands)
        end
    end
    unique_n_components = Set{Int}([])
    for (tf, estimands) ∈ tf_estimands
        # Number of generated estimands
        n_traits = 4
        n_treat_comb_order_2 = 5
        n_treat_comb_order_3 = 2
        @test length(estimands) == n_traits*(n_treat_comb_order_2+n_treat_comb_order_3)
        for Ψ ∈ estimands
            push!(unique_n_components, length(Ψ.args))
        end
    end
    # positivity constraint
    unique_n_components == Set([3, 9])
end

@testset "Test allele-independent: with positivity constraint" begin
    tmpdir = mktempdir()
    parsed_args = Dict(
        "allele-independent" => Dict{String, Any}(
            "config" => joinpath(TESTDIR, "data", "interaction_config.yaml"), 
            ),
        "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
        "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
        "call-threshold" => 0.8,  
        "%COMMAND%" => "allele-independent", 
        "bgen-prefix" => joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => joinpath(tmpdir, "final"), 
        "batch-size" => nothing,
        "positivity-constraint" => 0.01,
    )
    tmle_inputs(parsed_args)
    # Check dataset
    trait_data = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test sort(names(trait_data)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2", "RSID_2", "RSID_102", 
        "RSID_17", "RSID_198", "RSID_99"])
    @test size(trait_data) == (490, 16)
    # Check estimands
    tf_estimands = Dict(:TF1 => [], :TF2 => [])
    for file in readdir(tmpdir)
        if startswith(file, "final.TF1")
            append!(tf_estimands[:TF1], deserialize(joinpath(tmpdir, file)).estimands)
        elseif startswith(file, "final.TF2")
            append!(tf_estimands[:TF2], deserialize(joinpath(tmpdir, file)).estimands)
        end
    end
    unique_n_components = Set{Int}([])
    for (tf, estimands) ∈ tf_estimands
        # Number of generated estimands
        length(estimands) < 28
        for Ψ ∈ estimands
            push!(unique_n_components, length(Ψ.args))
        end
    end
    # positivity constraint
    unique_n_components == Set([1, 3, 5])
end

end

true