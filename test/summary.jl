module TestSummary

using Test
using TargeneCore
using CSV
using DataFrames

include("helper_fns.jl")


sieve_filename(prefix) = string(prefix, "_sieve_variance.hdf5")
summary_filename(prefix) = string(prefix, "_summary.csv")

function build_sieve_file(prefix)
    sieve_file = jldopen(sieve_filename(prefix), "w")
    sieve_file["STDERRORS"] = [0.3, nothing]
    sieve_file["SOURCEFILE_REPORTID_PAIRS"] = 
        ["$(prefix)_batch_1_Real.hdf5" => "2_1", "$(prefix)_batch_1_Real.hdf5" => "2_2"]
    close(sieve_file)
end

@testset "Test sieve_result_" begin
    @test all(x===missing for x in TargeneCore.sieve_results_(0.3, nothing))
    @test all(x isa Float64 for x in TargeneCore.sieve_results_(0.3, 0.1))
end

@testset "Test build_summary: with influence curve" begin
    grm_ids = TargeneCore.GRMIDs("data/grm/test.grm.id")
    prefix = "rs12_rs45"
    summaryfilename = summary_filename(prefix)
    build_results_files(grm_ids, prefix)
    # With sieve is false
    parsed_args = Dict(
        "prefix" => prefix,
        "out" => summary_filename(prefix),
        "sieve" => false
    )
    build_summary(parsed_args)

    summary_no_sieve = CSV.read(summaryfilename, DataFrame)

    rm(summaryfilename)

    # With sieve is true
    build_sieve_file(prefix)
    parsed_args = Dict(
        "prefix" => prefix,
        "out" => summary_filename(prefix),
        "sieve" => true
    )

    build_summary(parsed_args)
    
    summary = CSV.read(summaryfilename, DataFrame)

    @test summary.PHENOTYPE == ["cancer", "cancer", "height", "height", "bmi", "bmi"]
    @test summary.PHENOTYPE_TYPE == ["BINARY", "BINARY", "CONTINUOUS", "CONTINUOUS", "CONTINUOUS", "CONTINUOUS"]
    @test summary.QUERYNAME == ["QUERY_1", "QUERY_2", "QUERY_1", "QUERY_2", "QUERY_1", "QUERY_2"]
    @test summary.FILENAME_ORIGIN == ["rs12_rs45_batch_1_Bool.hdf5", 
                                      "rs12_rs45_batch_1_Bool.hdf5", 
                                      "rs12_rs45_batch_1_Real.hdf5", 
                                      "rs12_rs45_batch_1_Real.hdf5", 
                                      "rs12_rs45_batch_1_Real.hdf5", 
                                      "rs12_rs45_batch_1_Real.hdf5"]
    @test summary.REPORT_KEY == ["1_1", "1_2", "1_1", "1_2", "2_1", "2_2"]
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

    # Check the summary with and without sieve match on the same columns
    sievecolnames = names(TargeneCore.init_sieve("toto"))
    @test summary_no_sieve == summary[!, Not(sievecolnames)]

end

@testset "Test build_summary: without influence curve" begin
    grm_ids = TargeneCore.GRMIDs("data/grm/test.grm.id")
    prefix = "rs12_rs45"
    summaryfilename = summary_filename(prefix)
    build_results_files(grm_ids, prefix;save_ic=false)
    # With sieve is false
    parsed_args = Dict(
        "prefix" => prefix,
        "out" => summary_filename(prefix),
        "sieve" => false
    )
    build_summary(parsed_args)

    summary_no_sieve = CSV.read(summaryfilename, DataFrame)
    # Check the influence has indeed not been saved
    jldopen("rs12_rs45_batch_1_Bool.hdf5") do io
        @test io["TMLEREPORTS"]["1_1"].stderror isa Real
        @test io["TMLEREPORTS"]["1_1"].pvalue isa Real
        @test io["TMLEREPORTS"]["1_1"].confint isa Tuple
    end
    @test names(summary_no_sieve) == ["PHENOTYPE",
                                    "PHENOTYPE_TYPE",
                                    "QUERYNAME",
                                    "FILENAME_ORIGIN",
                                    "REPORT_KEY",
                                    "INITIAL_ESTIMATE",
                                    "ESTIMATE",
                                    "STDERR",
                                    "PVAL",
                                    "LWB",
                                    "UPB"]
    @test size(summary_no_sieve) == (6, 11)
    rm(summaryfilename)
    clean_hdf5estimates_files(prefix)
end


end

true