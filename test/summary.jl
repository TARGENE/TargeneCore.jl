module TestSummary

using Test
using TMLEEpistasis
using CSV
using DataFrames

include("helper_fns.jl")


sieve_filename(prefix) = string(prefix, "_sieve_variance.hdf5")
summary_filename(prefix) = string(prefix, "_summary.csv")

function build_sieve_file(prefix)
    sieve_file = jldopen(sieve_filename(prefix), "w")
    sieve_file["STDERRORS"] = [0.3, nothing]
    sieve_file["SOURCEFILE_REPORTID_PAIRS"] = 
        ["$(prefix)_batch_1_Real.hdf5" => 3, "$(prefix)_batch_1_Real.hdf5" => 4]
    close(sieve_file)
end

@testset "Test sieve_result" begin
    @test TMLEEpistasis.sieve_results(0.3, nothing) === 
        (missing, missing, missing, missing)
    
    @test all(x isa Float64 for x in TMLEEpistasis.sieve_results(0.3, 0.1))
end

@testset "Test build_summary" begin
    grm_ids = TMLEEpistasis.GRMIDs("data/grm/test.grm.id")
    prefix = "rs12_rs45"
    build_results_files(grm_ids, prefix; mode="QUERYREPORTS")
    build_sieve_file(prefix)
    parsed_args = Dict(
        "prefix" => prefix,
        "out" => summary_filename(prefix)
    )

    build_summary(parsed_args)
    
    summaryfilename = summary_filename(prefix)
    summary = CSV.read(summaryfilename, DataFrame)

    @test summary.PHENOTYPE == ["cancer", "cancer", "height", "height", "bmi", "bmi"]
    @test summary.QUERYNAME == ["QUERY_1", "QUERY_2", "QUERY_1", "QUERY_2", "QUERY_1", "QUERY_2"]
    @test eltype(summary.INITIAL_ESTIMATE) == Float64
    @test eltype(summary.ESTIMATE) == Float64
    @test eltype(summary.PVAL) == Float64
    @test eltype(summary.LWB) == Float64
    @test eltype(summary.UPB) == Float64
    @test summary.SIEVE_STDERR[5] == 0.3
    @test summary.SIEVE_PVAL[5] isa Float64
    @test summary.SIEVE_LWB[5] isa Float64
    @test summary.SIEVE_UPB[5] isa Float64

    clean_hdf5estimates_files(prefix)
    rm(sieve_filename(prefix))
    rm(summaryfilename)

end

end

true