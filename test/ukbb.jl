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

function test_base_serialization(queryreports)
    n_expected = 477
    queryreport = queryreports[1]
    @test queryreport isa TMLE.QueryReport
    test_queries((queryreport.query,),
        (Query(name="QUERY_1", case=(RSID_10 = "AG", RSID_100 = "AG"), control=(RSID_10 = "GG", RSID_100 = "GG")),)
    )
    @test size(queryreport.influence_curve, 1) == n_expected
    @test queryreport.estimate isa Real
    @test queryreport.initial_estimate isa Real

    # Check second queryreport
    queryreport = queryreports[2]
    @test queryreport isa TMLE.QueryReport
    test_queries((queryreport.query,),
        (Query(name="QUERY_2", case=(RSID_10 = "AG", RSID_100 = "AA"), control=(RSID_10 = "GG", RSID_100 = "GG")),)
    )
    @test size(queryreport.influence_curve, 1) == n_expected
    @test queryreport.estimate isa Real
    @test queryreport.initial_estimate isa Real
end

@testset "Test read_bgen function" begin
    # This file as an associated .bgi
    bgen_file = BGEN.datadir("example.8bits.bgen")
    b = TMLEEpistasis.read_bgen(bgen_file)
    @test b.idx isa BGEN.Index

    # This file as no associated .bgi
    bgen_file = BGEN.datadir("example.9bits.bgen")
    b = TMLEEpistasis.read_bgen(bgen_file)
    @test b.idx === nothing
end

@testset "Test samples_genotype" begin
    probabilities = [0.3 0.2 0.9;
                     0.5 0.2 0.05;
                     0.2 0.6 0.05]
    variant_genotypes = ["AA", "AT", "TT"]

    threshold = 0.9
    genotypes = TMLEEpistasis.samples_genotype(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2] === missing
    @test genotypes[3] == "AA"

    threshold = 0.55
    genotypes = TMLEEpistasis.samples_genotype(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1]  === missing
    @test genotypes[2] == "TT"
    @test genotypes[3] == "AA"

end

@testset "Test UKBBGenotypes function" begin
    build_query_file(threshold=0.95)
    queries = TMLEEpistasis.parse_queries(queryfile)
    genotypes = TMLEEpistasis.UKBBGenotypes(queryfile, queries)
    # I only look at the first 10 rows
    # SAMPLE_ID    
    @test genotypes[1:9, "SAMPLE_ID"] == ["sample_00$i" for i in 1:9]
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == repeat(["AG"], 10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== ["AG", "AA", "AG", missing, "AG", "AG", missing, "AG", "GG", "AG"])
    # Test column order
    @test DataFrames.names(genotypes) == ["SAMPLE_ID", "RSID_10", "RSID_100"]

    rm(queryfile)
end

@testset "Test variant_genotypes" begin
    b = TMLEEpistasis.read_bgen(BGEN.datadir("example.8bits.bgen"))
    queries = [
        Query(case=(RSID_10="AA", RSID_100="AA"), control=(RSID_10="GG", RSID_100="GG")),
        Query(case=(RSID_10="AA", RSID_100="AA"), control=(RSID_10="GG", RSID_100="GA"))
    ]
    # Defaults to "AG"
    v = variant_by_rsid(b, "RSID_10")
    @test TMLEEpistasis.variant_genotypes(v, queries) == ["AA", "AG", "GG"]
    # "GA" is specified by user so returned
    v = variant_by_rsid(b, "RSID_100")
    @test TMLEEpistasis.variant_genotypes(v, queries) == ["AA", "GA", "GG"]

end


@testset "Test preprocess" begin
    confounders = CSV.File(confoundersfile) |> DataFrame
    n = size(confounders)[1]
    genotypes = DataFrame(
        SAMPLE_ID = confounders.SAMPLE_ID, 
        RSID_1    = vcat(Vector{Missing}(missing, 10), collect(1:n-10)),
        RSID_2    = vcat(collect(1:n-10), Vector{Missing}(missing, 10))
        )

    # Continuous phenotypes
    phenotypes = CSV.File(phenotypefile, select=["eid", "continuous_phenotype"]) |> DataFrame
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes;verbosity=0)

    @test y == phenotypes.continuous_phenotype[11:n-10]

    n_expected = 480
    @test size(T) == (480, 2)
    @test size(W) == (480, 2)
    @test size(y, 1) == n_expected
    @test size(sample_ids, 1) == n_expected
    @test all(startswith(x, "sample") for x in sample_ids)

    @test DataFrames.names(T) == ["RSID_1", "RSID_2"]
    @test DataFrames.names(W) == ["PC1", "PC2"]

    @test T.RSID_1 == collect(1:n-20)
    @test T.RSID_1 isa CategoricalArray
    @test T.RSID_2 == collect(11:n-10)
    @test T.RSID_2 isa CategoricalArray

    # Binary phenotypes
    phenotypes = CSV.File(phenotypefile, select=["eid", "categorical_phenotype"]) |> DataFrame
    T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes;
                verbosity=0)

    @test y isa CategoricalArray{Bool,1,UInt32}
end


@testset "Test PhenotypeTMLEEpistasis with continuous target" begin
    # Only one continuous phenotype / machines not saved / no adaptive cv
    estimatorfile = joinpath("config", "tmle_config.toml")
    build_query_file()
    parsed_args = Dict(
        "phenotypes" => phenotypefile,
        "confounders" => confoundersfile,
        "queries" => queryfile,
        "estimator" => estimatorfile,
        "phenotypes-list" => phenotypelist_file_1,
        "verbosity" => 0,
        "adaptive-cv" => false,
        "save-full" => false
    )

    UKBBVariantRun(parsed_args)

    # Essential results
    hdf5_file = "RSID_10_RSID_100.hdf5"
    file = jldopen(hdf5_file)
    n_expected = 477
    @test size(file["continuous_phenotype"]["sample_ids"], 1) == n_expected
    test_base_serialization(file["continuous_phenotype"]["queryreports"])
    close(file)

    # Clean
    rm(hdf5_file)
    rm(parsed_args["queries"])
end


@testset "Test PhenotypeTMLEEpistasis with all phenotypes" begin
    estimatorfile = joinpath("config", "tmle_config.toml")
    build_query_file()
    parsed_args = Dict(
        "phenotypes" => phenotypefile,
        "confounders" => confoundersfile,
        "queries" => queryfile,
        "estimator" => estimatorfile,
        "phenotypes-list" => nothing,
        "verbosity" => 0,
        "adaptive-cv" => true,
        "save-full" => true
    )

    UKBBVariantRun(parsed_args)
    
    # Essential results
    hdf5_file = "RSID_10_RSID_100.hdf5"
    file = jldopen(hdf5_file)
    n_expected = 477
    @test size(file["continuous_phenotype"]["sample_ids"], 1) == n_expected
    @test size(file["categorical_phenotype"]["sample_ids"], 1) == n_expected
    test_base_serialization(file["continuous_phenotype"]["queryreports"])
    test_base_serialization(file["categorical_phenotype"]["queryreports"])
    close(file)

    # Machine
    jls_file = "RSID_10_RSID_100.jls"
    mach_file = open(jls_file)
    # First phenotype
    phenotype, tmle_mach = deserialize(mach_file)
    @test phenotype == :categorical_phenotype
    @test length(getqueryreports(tmle_mach)) == 2
    @test length(report(tmle_mach).G.cv_report) == 3
    @test length(report(tmle_mach).Q̅.cv_report) == 4
    test_data_has_been_removed(tmle_mach)
    # Second phenotype
    phenotype, tmle_mach = deserialize(mach_file)
    @test phenotype == :continuous_phenotype
    @test length(getqueryreports(tmle_mach)) == 2
    @test length(report(tmle_mach).G.cv_report) == 3
    @test length(report(tmle_mach).Q̅.cv_report) == 5
    test_data_has_been_removed(tmle_mach)
    # Adaptive CV 
    @test size(report(tmle_mach).Q̅.cv_report.HALRegressor_1.per_fold[1], 1) == 20

    # Only the continuous phenotype has been processed
    @test eof(mach_file)

    close(mach_file)

    # Clean
    rm(hdf5_file)
    rm(jls_file)
    rm(parsed_args["queries"])
end

end;

true