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

@testset "Test frequency_criterion" begin
    # Some combination is not present in the data so 
    # the criterion should return false whatever the threshold
    data = DataFrame(
        RSID_1=["AC", "AC", "AA", "AA", "AA"],
        RSID_2=["AC", "AA", "AA", "AA", "AA"],
        TRAIT=[1, 1, 0, 1, 1]
    )
    threshold = 0.
    @test TMLEEpistasis.frequency_criterion(data, 
        [:RSID_1, :RSID_2]; threshold=threshold) == false
    # Now with the same threshold it is true
    data = DataFrame(
        RSID_1=["AC", "AC", "AA", "AA", "AA", "AA"],
        RSID_2=["AC", "AA", "AA", "AA", "AA", "AC"],
        TRAIT=[1, 1, 0, 1, 1, 1]
    )
    @test TMLEEpistasis.frequency_criterion(data, 
        [:RSID_1, :RSID_2]; threshold=threshold) == true
    # Let's modify the threshold, lowest freq is 0.166
    threshold = 0.17
    @test TMLEEpistasis.frequency_criterion(data, 
        [:RSID_1, :RSID_2]; threshold=threshold) == false
    # including some phenotype, not all combinations are present
    threshold =0.
    @test TMLEEpistasis.frequency_criterion(data, 
        [:RSID_1, :RSID_2, :TRAIT]; threshold=threshold) == false

end

@testset "Test preprocess" begin
    n = 10
    base_sample_ids = 1:n

    confounders = DataFrame(rand(n, 3), :auto)
    confounders.SAMPLE_ID = base_sample_ids

    genotypes = DataFrame(
        SAMPLE_ID = base_sample_ids, 
        RSID_1    = ["AC", "AC", "AC", "AC", "AC", "CC", "CC", "CC", "CC", "CC"],
        RSID_2    = ["AG", "AA", "AA", "AA", "AA", "AG", "AG", "AA", "AA", "AA"]
        )

    queries = [TMLE.Query(case=(RSID_2="AA", RSID_1="CC"), control=(RSID_2="AG", RSID_1="AC"))]
    phenotypes = DataFrame(
        SAMPLE_ID = base_sample_ids,
        Y1 = rand(n),
        Y2 = rand(n)
    )
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, Real, queries, freq_threhsold=0.01)
    
    # Check the order of columns has been reversed
    @test T == genotypes[!, [:RSID_2, :RSID_1]]
    @test W == confounders[!, Not("SAMPLE_ID")]
    @test y == phenotypes[!, Not("SAMPLE_ID")]
    @test sample_ids == 
            Dict(
                "Y1" => string.(base_sample_ids),
                "Y2" => string.(base_sample_ids),
            )
    
    # frequency threshold not passed by genotype due to missing combination
    genotypes = DataFrame(
        SAMPLE_ID = base_sample_ids, 
        RSID_1    = ["AC", "AC", "AC", "AC", "AC", "CC", "CC", "CC", "CC", "CC"],
        RSID_2    = [missing, "AA", "AA", "AA", "AA", "AG", "AG", "AA", "AA", "AA"]
        )
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, Real, queries, freq_threhsold=0.01)
    @test all((T, W, y, sample_ids) .=== (nothing, nothing, nothing, nothing))

    # frequency threshold not passed by genotype due to poor frequency
    genotypes = DataFrame(
        SAMPLE_ID = base_sample_ids, 
        RSID_1    = ["AC", "AC", "AC", "AC", "AC", "CC", "CC", "CC", "CC", "CC"],
        RSID_2    = ["AG", "AA", "AA", "AA", "AA", "AG", "AG", "AA", "AA", "AA"]
        )
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, Real, queries, freq_threhsold=0.2)
    @test all((T, W, y, sample_ids) .=== (nothing, nothing, nothing, nothing))

    # frequency threshold not passed by some phenotype
    genotypes = DataFrame(
        SAMPLE_ID = base_sample_ids, 
        RSID_1    = ["AC", "AC", "AC", "AC", "AC", "CC", "CC", "CC", "CC", "CC"],
        RSID_2    = ["AG", "AG", "AA", "AA", "AA", "AG", "AG", "AA", "AA", "AA"]
        )
    phenotypes = DataFrame(
        SAMPLE_ID = base_sample_ids,
        Y1 = [0, 1, 0, 1, 0, 1, 0, 0, 1, 0],
        Y2 = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0]
    )
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, Bool, queries, freq_threhsold=0.0)
    @test T == genotypes[!, [:RSID_2, :RSID_1]]
    @test W == confounders[!, Not("SAMPLE_ID")]
    @test y == phenotypes[!, [:Y1]]
    @test y[!, :Y1] isa CategoricalVector
    @test sample_ids == 
            Dict(
                "Y1" => string.(base_sample_ids),
            )
    # frequency threshold not passed by all phenotype
    phenotypes = DataFrame(
        SAMPLE_ID = base_sample_ids,
        Y2 = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0]
    )
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, Bool, queries, freq_threhsold=0.0)
    @test all((T, W, y, sample_ids) .=== (nothing, nothing, nothing, nothing))

    # Checking joining on sample_ids
    genotypes = DataFrame(
        SAMPLE_ID = base_sample_ids, 
        RSID_1    = ["AC", "AC", missing, "AC", "AC", "CC", "CC", "CC", "CC", "CC"],
        RSID_2    = ["AG", "AG", "AA", "AA", "AA", "AG", "AG", "AA", "AA", "AA"]
        )
    phenotypes = DataFrame(
        SAMPLE_ID = base_sample_ids,
        Y1 = [0, 1, 0, 1, 0, 1, 0, 0, 1, missing],
    )
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, Bool, queries, freq_threhsold=0.0)

    @test size(T) == (9, 2)
    # For y the missing value hasn't been dropped yet, it will be dropped during TMLE
    @test size(y) == (9, 1)
    @test size(W) == (9, 3)
    @test sample_ids == Dict("Y1" => ["1", "2", "4", "5", "6", "7", "8", "9"])

end


@testset "Test UKBBVariantRun with continuous target" begin
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
        "target-type" => "Real",
        "min-freq" => 0.01
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


@testset "Test UKBBVariantRun with binary targets" begin
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
        "target-type" => "Bool",
        "min-freq" => 0.01
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

@testset "Test UKBBVariantRun with high freq threhsold" begin
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
        "target-type" => "Bool",
        "min-freq" => 0.05
    )

    UKBBVariantRun(parsed_args)
    
    # Essential results
    jldopen(parsed_args["out"]) do io
        @test isempty(keys(io["TMLEREPORTS"]))
    end
    # Clean
    rm(parsed_args["out"])
    rm(parsed_args["queries"])
end

end;

true