module TestUtils

using Test
using TMLEEpistasis
using MLJ
using TOML
using BGEN
using CSV
using DataFrames
using Serialization
using StableRNGs
using TMLE

include("helper_fns.jl")

LinearRegressor = @load LinearRegressor pkg=MLJLinearModels verbosity=0
LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity=0


@testset "Test parse_queries" begin
    build_query_file()
    queries = TMLEEpistasis.parse_queries(queryfile)
    @test queries == [
        "QUERY_1" => (RSID_10 = ["AG", "GG"], RSID_100 = ["AG", "GG"]),
        "QUERY_2" => (RSID_10 = ["AG", "GG"], RSID_100 = ["AA", "GG"])
    ]
    rm(queryfile)

    @test TMLEEpistasis.querystring(queries[1][2]) == "RSID_10: AG -> GG & RSID_100: AG -> GG"
end

@testset "Test phenotypes parsing" begin
    allnames = TMLEEpistasis.phenotypesnames(phenotypefile)
    @test allnames == [:categorical_phenotype, :continuous_phenotype]
end

@testset "Test set_cv_folds!" begin
    rng = StableRNG(123)
    tmle_config = joinpath("config", "tmle_config.toml")
    build_query_file()
    queries = TMLEEpistasis.parse_queries(queryfile)
    tmles =  TMLEEpistasis.estimators_from_toml(TOML.parsefile(tmle_config), queries)

    # Continuous y
    tmle = tmles["continuous"]
    @test tmle.Q̅.resampling.nfolds == 2

    ## adaptive_cv = false
    y = rand(100)
    TMLEEpistasis.set_cv_folds!(tmle, y, learner=:Q̅, adaptive_cv=false)
    @test tmle.Q̅.resampling.nfolds == 2
    ## neff = n = 100 => 20 folds
    TMLEEpistasis.set_cv_folds!(tmle, y, learner=:Q̅, adaptive_cv=true)
    @test tmle.Q̅.resampling.nfolds == 20
    ## neff = n = 20_000 => 3 folds
    y = rand(20_000)
    TMLEEpistasis.set_cv_folds!(tmle, y, learner=:Q̅, adaptive_cv=true)
    @test tmle.Q̅.resampling.nfolds == 3


    # Categorical y
    tmle = tmles["binary"]
    @test tmle.Q̅.resampling.nfolds == 2
    y = categorical(["a", "a", "a", "b", "b", "c", "c"])
    ## neff < 30 => nfolds = 5*neff = 7
    TMLEEpistasis.set_cv_folds!(tmle, y, learner=:Q̅, adaptive_cv=true)
    @test tmle.Q̅.resampling.nfolds == 7
    ## neff = 2500 => 10
    y = categorical(vcat(repeat([true], 50_000), repeat([false], 500)))
    TMLEEpistasis.set_cv_folds!(tmle, y, learner=:Q̅, adaptive_cv=true)
    @test tmle.Q̅.resampling.nfolds == 10

    rm(queryfile)

end

@testset "Test writeresults" begin
    # Build a TMLEstimator
    query = (t₁=["CG", "GG"], t₂=["TT", "TA"])
    Q = LinearRegressor()
    G = FullCategoricalJoint(LogisticClassifier())
    tmle = TMLEstimator(Q, G, query)
    # Fit with arbitrary data 
    n = 100
    W = (x₁=rand(n), )
    T = (
        t₁=categorical(rand(["CG", "GG"], n)), 
        t₂=categorical(rand(["TT", "TA"], n))
    )
    y = rand(n)
    mach = machine(tmle, T, W, y)
    fit!(mach, verbosity=0)
    # Save the results
    outfile = "results.bin"
    open(outfile, "w") do file
        TMLEEpistasis.writeresults(file, mach, :cancer, full=false)
        TMLEEpistasis.writeresults(file, mach, "hair_color", full=true)
    end
    # Load back
    open(outfile) do file
        phenotype, reports = deserialize(file)
        @test phenotype == :cancer
        @test all(qr isa TMLE.QueryReport for qr in reports)

        phenotype, mach = deserialize(file)
        @test phenotype == "hair_color"
        @test mach isa Machine
    end

    rm(outfile)
end

end;

true