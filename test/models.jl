module TestModels

using Test
using TMLEEpistasis
using MLJBase
using CategoricalDistributions


@testset "Test InteractionTransformer" begin
    X = (rs1234=[1, 2, 3], rs455=[4, 5, 6], rs4489=[7, 8, 9], rstoto=[1, 2, 3])
    t = TMLEEpistasis.InteractionTransformer(r"^rs[0-9]+")
    mach = machine(t, X)
    fit!(mach, verbosity=0)
    Xt = MLJBase.transform(mach, X)

    @test Xt == (
        rs1234 = [1, 2, 3],
        rs455 = [4, 5, 6],
        rs4489 = [7, 8, 9],
        rstoto = [1, 2, 3],
        rs1234_rs455 = [4.0, 10.0, 18.0],
        rs1234_rs4489 = [7.0, 16.0, 27.0],
        rs455_rs4489 = [28.0, 40.0, 54.0]
    )
    @test mach.fitresult.ninter == 3
    @test mach.fitresult.interaction_pairs == [:rs1234 => :rs455, :rs1234 => :rs4489, :rs455 => :rs4489]

end

@testset "Test InteractionLM" begin
    n = 100
    X = (rs1234=rand(n), rs455=rand(n), rstoto=rand(n))
    # Classifier
    y = categorical(rand([0,1], n))
    mach = machine(
        TMLEEpistasis.InteractionLMClassifier(),
        X,
        y
    )
    fit!(mach, verbosity=0)
    fp = fitted_params(mach)
    @test fp.interaction_transformer.fitresult.interaction_pairs == [:rs1234 => :rs455]
    @test predict(mach) isa UnivariateFiniteArray{Multiclass{2}, Int64, UInt32, Float64, 1}
    # Regressor 
    y = rand(n)
    mach = machine(
        TMLEEpistasis.InteractionLMRegressor(),
        X,
        y
    )
    fit!(mach, verbosity=0)
    fp = fitted_params(mach)
    @test fp.interaction_transformer.fitresult.interaction_pairs == [:rs1234 => :rs455]
    @test predict(mach) isa Vector{Float64}

    @test target_scitype(TMLEEpistasis.InteractionLMRegressor()) == AbstractVector{<:MLJBase.Continuous}
    @test target_scitype(TMLEEpistasis.InteractionLMClassifier()) == AbstractVector{<:MLJBase.Finite}

end

end
true