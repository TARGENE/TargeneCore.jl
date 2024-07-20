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

end

true