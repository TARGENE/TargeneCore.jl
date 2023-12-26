module TestAlleleIndependeEstimands

using Test
using TargeneCore
using Arrow
using DataFrames
using Serialization

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "tmle_inputs", "test_utils.jl"))

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
    for (tf, estimands) ∈ tf_estimands
        # Number of generated estimands
        ntraits, nbQTLs, neQTLs = 4, 2, 1
        @test length(estimands) == ntraits*nbQTLs*neQTLs
        for Ψ ∈ estimands
            # No positivity constraint
            @test length(Ψ.args) == 9
        end
    end
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
    for (tf, estimands) ∈ tf_estimands
        # Number of generated estimands
        @test length(estimands) == 4 < 8
        for Ψ ∈ estimands
            # No positivity constraint
            @test length(Ψ.args) == 1 < 9
        end
    end
end

end

true