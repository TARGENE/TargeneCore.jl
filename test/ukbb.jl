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
    n_expected = 466
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
    # Let's remove the first 10 lines of the genotypes
    genotypes = DataFrame(
        SAMPLE_ID = confounders.SAMPLE_ID, 
        RSID_1    = collect(1:n),
        RSID_2    = vcat(collect(1:10), missing, collect(12:n))
        )[11:end, :]

    # Continuous phenotypes
    for (phen, target_type) in (("continuous_phenotype", Real), ("categorical_phenotype", Bool))
        phenotypes = CSV.File(phenotypefile, select=["SAMPLE_ID", phen]) |> DataFrame
        T, W, y, sample_ids = TMLEEpistasis.preprocess(genotypes, confounders, phenotypes, target_type)
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
        @test length(sample_ids[phen]) == 478
        if phen == "continuous_phenotype"
            @test y[1:end-1, :] == DataFrame(phen => phenotypes[12:489, phen])
            y[end, 1] === missing
            # sample 490 is missing
            @test sample_ids[phen][end] == "sample_489"
        else
            @test y[1:end-3, :] == DataFrame(phen => categorical(phenotypes[12:487, phen]))
            y[end-2, 1] === missing
            # sample 489 is missing
            @test sample_ids[phen][end-1] == "sample_488"
        end
    end
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
        "save-full" => false,
        "target-type" => "Real"
    )

    UKBBVariantRun(parsed_args)
    # Essential results
    outfilename = "RSID_10_RSID_100.hdf5"
    file = jldopen(outfilename)
    n_expected = 466
    @test size(file["SAMPLE_IDS"]["continuous_phenotype"], 1) == n_expected
    test_base_serialization(file["QUERYREPORTS"])
    close(file)

    # Clean
    rm(outfilename)
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
        "phenotypes-list" => phenotypelist_file_2,
        "verbosity" => 0,
        "adaptive-cv" => true,
        "save-full" => true,
        "target-type" => "Bool"
    )

    UKBBVariantRun(parsed_args)
    
    # Essential results
    outfilename = "RSID_10_RSID_100.hdf5"
    file = jldopen(outfilename)
    n_expected = 466
    @test size(file["SAMPLE_IDS"]["categorical_phenotype"], 1) == n_expected
    mach = file["MACHINE"]
    test_base_serialization(queryreports(mach))

    @test length(report(mach).G.cv_report) == 3
    @test length(report(mach).Q̅.cv_report) == 4

    # Adaptive CV 
    @test size(report(mach).Q̅.cv_report.HALClassifier_1.per_fold[1], 1) == 20

    # Clean
    rm(outfilename)
    rm(parsed_args["queries"])
end

end;

true