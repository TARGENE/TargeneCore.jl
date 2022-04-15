module TestsUKBB

using Test
using BGEN
using CSV
using DataFrames
using TMLEEpistasis
using CategoricalArrays
using TOML
using Serialization
using TMLE
using JLD2
using MLJBase

include("helper_fns.jl")

function test_base_serialization(queryreports, n_expected)
    queryreport = queryreports[1]
    @test queryreport isa TMLE.Report
    test_queries((queryreport.query,),
        (Query(name="QUERY_1", case=(RSID_10 = "AG", RSID_100 = "AG"), control=(RSID_10 = "GG", RSID_100 = "GG")),)
    )
    @test size(queryreport.influence_curve, 1) == n_expected
    @test queryreport.estimate isa Real
    @test queryreport.initial_estimate isa Real

    # Check second queryreport
    queryreport = queryreports[2]
    @test queryreport isa TMLE.Report
    test_queries((queryreport.query,),
        (Query(name="QUERY_2", case=(RSID_10 = "AG", RSID_100 = "AA"), control=(RSID_10 = "GG", RSID_100 = "GG")),)
    )
    @test size(queryreport.influence_curve, 1) == n_expected
    @test queryreport.estimate isa Real
    @test queryreport.initial_estimate isa Real
end

@testset "Test preprocess" begin
    confounders = CSV.File(confoundersfile) |> DataFrame
    n = size(confounders)[1]
    # Let's remove the first 10 lines of the genotypes
    genotypes = DataFrame(
        SAMPLE_ID = confounders.SAMPLE_ID, 
        RSID_1    = collect(1:n),
        RSID_2    = vcat(collect(1:10), missing, collect(12:n))
        )[11:end, :]

    # Continuous phenotypes
    phenotypes = CSV.File(continuous_phenotypefile) |> DataFrame
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, Real)
    @test T == DataFrame(
                RSID_1 = categorical(genotypes[2:end-10, "RSID_1"]),
                RSID_2 = categorical(genotypes[2:end-10, "RSID_2"])
        )
    @test typeof(T.RSID_2) == CategoricalArray{
        Int64, 
        1, 
        UInt32, 
        Int64, 
        CategoricalValue{Int64, UInt32}, Union{}}
    @test W == confounders[12:end-10, ["PC1", "PC2"]]
    @test length(sample_ids["CONTINUOUS_1"]) == 478
    @test length(sample_ids["CONTINUOUS_2"]) == 477
    @test all(x === missing for x in y[478:479, "CONTINUOUS_2"])
    @test sample_ids["CONTINUOUS_2"][end] == "sample_488"  
    @test sample_ids["CONTINUOUS_1"][end-1] == "sample_488"  

    # Binary phenotypes
    phenotypes = CSV.File(binary_phenotypefile) |> DataFrame
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, Bool)
    @test T == DataFrame(
                RSID_1 = categorical(genotypes[2:end-10, "RSID_1"]),
                RSID_2 = categorical(genotypes[2:end-10, "RSID_2"])
        )
    @test typeof(T.RSID_2) == CategoricalArray{
        Int64, 
        1, 
        UInt32, 
        Int64, 
        CategoricalValue{Int64, UInt32}, Union{}}
    @test W == confounders[12:end-10, ["PC1", "PC2"]]
    @test length(sample_ids["BINARY_1"]) == 479
    @test length(sample_ids["BINARY_2"]) == 477
    @test all(x === missing for x in y[478:479, "BINARY_2"])
    @test sample_ids["BINARY_1"][end] == "sample_490"  
    @test sample_ids["BINARY_2"][end-1] == "sample_487"  
end


@testset "Test PhenotypeTMLEEpistasis with continuous target" begin
    # Only one continuous phenotype / machines not saved / no adaptive cv
    estimatorfile = joinpath("config", "tmle_config.toml")
    build_query_file()
    parsed_args = Dict(
        "phenotypes" => continuous_phenotypefile,
        "confounders" => confoundersfile,
        "queries" => queryfile,
        "estimator" => estimatorfile,
        "out" => "RSID_10_RSID_100.hdf5",
        "phenotypes-list" => phenotypelist_file,
        "verbosity" => 0,
        "adaptive-cv" => false,
        "save-full" => false,
        "target-type" => "Real"
    )

    UKBBVariantRun(parsed_args)
    # Essential results
    file = jldopen(parsed_args["out"])
    n_expected = 466
    @test size(file["SAMPLE_IDS"]["CONTINUOUS_1"], 1) == n_expected
    test_base_serialization(file["QUERYREPORTS"], n_expected)
    close(file)

    # Clean
    rm(parsed_args["out"])
    rm(parsed_args["queries"])
end


@testset "Test PhenotypeTMLEEpistasis with all phenotypes" begin
    estimatorfile = joinpath("config", "tmle_config.toml")
    build_query_file()
    parsed_args = Dict(
        "phenotypes" => binary_phenotypefile,
        "confounders" => confoundersfile,
        "queries" => queryfile,
        "estimator" => estimatorfile,
        "out" => "RSID_10_RSID_100.hdf5",
        "phenotypes-list" => nothing,
        "verbosity" => 0,
        "adaptive-cv" => true,
        "save-full" => true,
        "target-type" => "Bool"
    )

    UKBBVariantRun(parsed_args)
    
    # Essential results
    file = jldopen(parsed_args["out"])

    @test size(file["SAMPLE_IDS"]["BINARY_1"], 1) == 467
    @test size(file["SAMPLE_IDS"]["BINARY_2"], 1) == 465
    mach = file["MACHINE"]
    queryreports_ = queryreports(mach)
    test_base_serialization(filter(x -> x.target_name == :BINARY_1, queryreports_), 467)
    test_base_serialization(filter(x -> x.target_name == :BINARY_2, queryreports_), 465)
    @test length(report(mach).G.cv_report) == 3
    @test length(report(mach).Q̅[1].cv_report) == 4
    @test length(report(mach).Q̅[2].cv_report) == 4
    # Adaptive CV 
    @test size(report(mach).Q̅[1].cv_report.HALClassifier_1.per_fold[1], 1) == 20

    # Clean
    rm(parsed_args["out"])
    rm(parsed_args["queries"])
end

end;

true