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

include("helper_fns.jl")


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

    @test TMLEEpistasis.phenotypes_list(nothing, [:categorical_phenotype], allnames) == [:continuous_phenotype]
    @test TMLEEpistasis.phenotypes_list(nothing, [], allnames) == allnames

    @test TMLEEpistasis.phenotypes_list("data/phen_list_1.csv", [], allnames) == [:continuous_phenotype]
    @test TMLEEpistasis.phenotypes_list("data/phen_list_2.csv", [:continuous_phenotype], allnames) == [:categorical_phenotype]
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

end;

true