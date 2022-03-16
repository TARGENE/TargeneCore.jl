module TestUtils

using Test
using TMLEEpistasis
using MLJBase
using TOML
using BGEN
using DataFrames
using Serialization
using TMLE
using MLJLinearModels

include("helper_fns.jl")

@testset "Test parse_queries" begin
    build_query_file()
    queries = TMLEEpistasis.parse_queries(queryfile)
    expected_queries = [
        Query(case=(RSID_10="AG", RSID_100="AG"), control=(RSID_10="GG", RSID_100="GG"), name="QUERY_1"),
        Query(case=(RSID_10="AG", RSID_100="AA"), control=(RSID_10="GG", RSID_100="GG"), name="QUERY_2")
    ]
    test_queries(queries, expected_queries)

    rm(queryfile)
end

@testset "Test out filenames" begin
    query = Query(case=(RSID_10="AG", RSID_100="AG"), control=(RSID_10="GG", RSID_100="GG"), name="QUERY_1")
    @test TMLEEpistasis.hdf5filename(query) == "RSID_10_RSID_100.hdf5"
    @test TMLEEpistasis.jlsfilename(query) == "RSID_10_RSID_100.jls"

    query = Query(case=(RSID_10="AG",), control=(RSID_10="GG",), name="QUERY_1")
    @test TMLEEpistasis.hdf5filename(query) == "RSID_10.hdf5"
    @test TMLEEpistasis.jlsfilename(query) == "RSID_10.jls"
end

@testset "Test phenotypes parsing" begin
    # Fallback when no list is specified
    @test nothing === TMLEEpistasis.phenotypesnames(nothing)
    # Test with a list of phenotypes
    @test ["FID", "SAMPLE_ID", "continuous_phenotype"] == TMLEEpistasis.phenotypesnames(phenotypelist_file_1)
    @test ["FID", "SAMPLE_ID", "categorical_phenotype"] == TMLEEpistasis.phenotypesnames(phenotypelist_file_2)

    # phenotypes loading
    phenotypes = TMLEEpistasis.load_phenotypes(phenotypefile, nothing)
    @test size(phenotypes) == (490, 4)
    @test names(phenotypes) == ["FID", "SAMPLE_ID", "categorical_phenotype", "continuous_phenotype"]
    phenotypes = TMLEEpistasis.load_phenotypes(phenotypefile, phenotypelist_file_1)
    @test size(phenotypes) == (490, 3)
    @test names(phenotypes) == ["FID", "SAMPLE_ID", "continuous_phenotype"]
end

@testset "Test serializable!" begin
    tmle_config = joinpath("config", "tmle_config.toml")
    build_query_file()
    queries = TMLEEpistasis.parse_queries(queryfile)
    tmle =  TMLEEpistasis.estimator_from_toml(TOML.parsefile(tmle_config), queries, Bool; adaptive_cv=false)
    n = 100
    y = categorical(rand(Bool, n))
    T = (
        RSID_10=categorical(rand(["AA","AG", "GG"], n)),
        RSID_100=categorical(rand(["AA","AG", "GG"], n))
        )
    W = (w₁ = rand(n),)
    mach = machine(tmle, T, W, y)
    fit!(mach, verbosity=0)
    TMLEEpistasis.serializable!(mach)
    test_data_has_been_removed(mach)
    rm(queryfile)
end


@testset "Test AdaptiveCV" begin
    # Continuous target
    cv = TMLEEpistasis.AdaptiveCV(CV())
    n_sample_to_nfolds = ((10, 10), (200, 20), (1000, 10), (5000, 5), (20000, 3))
    for (n, expected_nfolds) in n_sample_to_nfolds
        y = rand(n)
        ttp = MLJBase.train_test_pairs(cv, 1:n, y)
        @test length(ttp) == expected_nfolds
        @test ttp == MLJBase.train_test_pairs(CV(nfolds=expected_nfolds), 1:n, y)
    end
    # Categorical target
    cv = TMLEEpistasis.AdaptiveCV(StratifiedCV(nfolds=2))
    y = categorical(["a", "a", "a", "b", "b", "c", "c"])
    @test TMLEEpistasis.countuniques(y) == [3, 2, 2]
    ## neff < 30 => nfolds = 5*neff = 7
    ttp = MLJBase.train_test_pairs(cv, 1:7, y)
    @test length(ttp) == 7
    @test ttp == MLJBase.train_test_pairs(StratifiedCV(nfolds=7), 1:7, y)
    
    ## neff = 2500 => 10
    n = 50_500
    y = categorical(vcat(repeat([true], 50_000), repeat([false], 500)))
    @test TMLEEpistasis.countuniques(y) == [50_000, 500]
    ttp = MLJBase.train_test_pairs(cv, 1:n, y)
    @test length(ttp)== 10
    @test ttp == MLJBase.train_test_pairs(StratifiedCV(nfolds=10), 1:n, y)

end

end;

true