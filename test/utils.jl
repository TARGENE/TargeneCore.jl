module TestUtils

using Test
using TargeneCore

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "testutils.jl"))

@testset "Test treatment_variables and outcome_variables" begin
    estimates = make_estimates()
    # Classic estimand
    Ψ = estimates[1].TMLE.estimand
    @test Ψ isa TMLE.StatisticalIATE
    @test TargeneCore.treatment_variables(Ψ) == [:RSID_103, :rs10043934]
    @test TargeneCore.outcome_variables(Ψ) == [Ψ.outcome]
    # Composed estimand
    Ψ = estimates[3].TMLE.estimand
    @test Ψ isa JointEstimand
    @test TargeneCore.treatment_variables(Ψ) == [:RSID_103, :rs10043934]
    @test TargeneCore.outcome_variables(Ψ) == [Symbol("High light scatter reticulocyte percentage")]
end

@testset "Test read_significant_results" begin
    prefix = "tmle_output"
    estimates = make_estimates()
    save(estimates;prefix=prefix)

    threshold = 1.
    for ext in ("jls", "json", "hdf5")
        results = TargeneCore.read_significant_results(string(prefix, "." , ext); threshold=threshold)
        length(results) == 4
        for index ∈ 1:4
            results[index] == estimates[index].TMLE.estimand
        end
    end

    threshold = 1e-9
    for ext in ("jls", "json", "hdf5")
        results = TargeneCore.read_significant_results(string(prefix, "." ,ext); threshold=threshold)
        @test length(results) == 2
    end

    threshold = 1e-100
    for ext in ("jls", "json", "hdf5")
        results = TargeneCore.read_significant_results(string(prefix, "." ,ext); threshold=threshold)
        @test length(results) == 0
    end

    clean()
end

end

true