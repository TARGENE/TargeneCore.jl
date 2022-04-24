module TestsUKBB

using Test
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

function test_base_serialization(tmle_reports, n_expected; phenotype_id=1)
    tmle_report = tmle_reports["$(phenotype_id)_1"]
    @test tmle_report isa TMLE.TMLEReport
    test_queries((tmle_report.query,),
        (Query(name="QUERY_1", case=(RSID_10 = "AG", RSID_100 = "AG"), control=(RSID_10 = "GG", RSID_100 = "GG")),)
    )
    @test size(tmle_report.influence_curve, 1) == n_expected
    @test tmle_report.estimate isa Real
    @test tmle_report.initial_estimate isa Real

    # Check second tmle_report
    tmle_report = tmle_reports["$(phenotype_id)_2"]
    @test tmle_report isa TMLE.TMLEReport
    test_queries((tmle_report.query,),
        (Query(name="QUERY_2", case=(RSID_10 = "AG", RSID_100 = "AA"), control=(RSID_10 = "GG", RSID_100 = "GG")),)
    )
    @test size(tmle_report.influence_curve, 1) == n_expected
    @test tmle_report.estimate isa Real
    @test tmle_report.initial_estimate isa Real
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

    queries = [TMLE.Query(case=(RSID_2=2, RSID_1=1), control=(RSID_2=2, RSID_1=2))]
    # Continuous phenotypes
    phenotypes = CSV.File(continuous_phenotypefile) |> DataFrame
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, Real, queries)
    @test T == DataFrame(
                RSID_2 = categorical(genotypes[2:end-10, "RSID_2"]),
                RSID_1 = categorical(genotypes[2:end-10, "RSID_1"])
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
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, Bool, queries)
    @test T == DataFrame(
        RSID_2 = categorical(genotypes[2:end-10, "RSID_2"]),
        RSID_1 = categorical(genotypes[2:end-10, "RSID_1"])
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
        "genotypes" => genotypesfile,
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
    n_expected = 488
    @test size(file["SAMPLE_IDS"]["CONTINUOUS_1"], 1) == n_expected
    test_base_serialization(file["TMLEREPORTS"], n_expected)
    close(file)

    # Clean
    rm(parsed_args["out"])
    rm(parsed_args["queries"])
end


@testset "Test PhenotypeTMLEEpistasis with all phenotypes" begin
    estimatorfile = joinpath("config", "tmle_config.toml")
    build_query_file()
    parsed_args = Dict(
        "genotypes" => genotypesfile,
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

    @test size(file["SAMPLE_IDS"]["BINARY_1"], 1) == 489
    @test size(file["SAMPLE_IDS"]["BINARY_2"], 1) == 487
    
    tmlereports = file["TMLEREPORTS"]
    test_base_serialization(tmlereports, 489, phenotype_id=1)
    test_base_serialization(tmlereports, 487, phenotype_id=2)

    machines = file["MACHINES"]
    Gmach = machines["G"]
    @test length(report(Gmach).cv_report) == 3
    Qmach₁ = machines["Q_1"]
    @test length(report(Qmach₁).cv_report) == 4
    Qmach₂ = machines["Q_2"]
    @test length(report(Qmach₂).cv_report) == 4
    # Adaptive CV 
    @test size(report(Qmach₁).cv_report.HALClassifier_1.per_fold[1], 1) == 20

    # Clean
    rm(parsed_args["out"])
    rm(parsed_args["queries"])
end

end;

true