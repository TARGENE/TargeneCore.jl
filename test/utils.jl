module TestUtils

using Test
using TMLEEpistasis
using MLJ
using TOML
using BGEN
using CSV
using DataFrames

include("helper_fns.jl")

@testset "Test retrieve coefficients" begin
    estimatorfile = joinpath("config", "tmle_config.toml")
    tmle_config = TOML.parsefile(estimatorfile)
    queries = ["query_1" => (t₁=[1, 0], t₂=[1, 0]),]
    tmle = TMLEEpistasis.estimators_from_toml(tmle_config, queries, TMLEEpistasis.PhenotypeTMLEEpistasis)

    n = 1000
    T = (t₁=categorical(rand([1, 0], n)), t₂=categorical(rand([1, 0], n)))
    W = (w₁=rand(n),)
    y = rand(n)
    mach = machine(tmle["continuous"], T, W, y)
    fit!(mach, verbosity=0)

    coefs = TMLEEpistasis.Qstack_coefs(mach)
    @test coefs isa Vector
    @test length(coefs) == 4

    repr = TMLEEpistasis.repr_Qstack_coefs(mach)
    @test repr isa String
    str_reprs = split(repr, ", ")
    modelstr = join([split(v, " => ")[1] for v in str_reprs], ",")
    @test modelstr == "HALRegressor,XGBoostRegressor,XGBoostRegressor,InteractionLMRegressor"

end

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

    # Test init_or_retrieve_results for PhenotypeTMLEEpistasis
    outfile = "temp_results.csv"
    @test TMLEEpistasis.init_or_retrieve_results(outfile, TMLEEpistasis.PhenotypeTMLEEpistasis) == Set()
    init_df = CSV.File(outfile) |> DataFrame
    @test names(init_df) == ["PHENOTYPE",
                            "QUERYNAME",
                            "QUERYSTRING",
                            "ESTIMATE",
                            "PVALUE",
                            "LOWER_BOUND",
                            "UPPER_BOUND",
                            "STD_ERROR",
                            "QSTACK_COEFS",
                            "NROWS",
                            "TCOUNTS"]

    CSV.write(outfile, DataFrame(PHENOTYPE=["toto", "tata", "toto"]))
    @test TMLEEpistasis.init_or_retrieve_results(outfile, TMLEEpistasis.PhenotypeTMLEEpistasis) == Set([:toto, :tata])
    rm(outfile)

    # Test init_or_retrieve_results for PhenotypeCrossValidation
    @test TMLEEpistasis.init_or_retrieve_results(outfile, TMLEEpistasis.PhenotypeCrossValidation) == Set()
    init_df = CSV.File(outfile) |> DataFrame
    @test names(init_df) == ["PHENOTYPE",
                            "VARIANTS",
                            "Q_RESULTSTRING",
                            "G_RESULTSTRING"]
    CSV.write(outfile, DataFrame(PHENOTYPE=["toto", "tata", "toto"]))
    @test TMLEEpistasis.init_or_retrieve_results(outfile, TMLEEpistasis.PhenotypeCrossValidation) == Set([:toto, :tata])
    rm(outfile)
end

@testset "Test TMLE REPORT PARSING" begin
    # Only one query, the reports is just a NamedTuple
    reports = (pvalue=0.05,)
    @test TMLEEpistasis.get_query_report(reports, 4) == reports
    reports = [(pvalue=0.05,), (pvalue=0.01,)]
    @test TMLEEpistasis.get_query_report(reports, 1) == reports[1]
    @test TMLEEpistasis.get_query_report(reports, 2) == reports[2]

end


@testset "Test encode" begin
    # Multiple treatments combined to joint treatment
    T = DataFrame(
        t₁=categorical(["AG", "GG", "AA", "AA", "AG"]),
        t₂=categorical(["CT", "CC", "TT", "CT", "CT"]),
        )

    t_target = TMLEEpistasis.encode(T)
    @test levels(t_target) == [1, 2, 9, 3, 5, 6, 8, 7, 4]
    @test t_target == categorical([5, 3, 7, 4, 5])

    # One treatment just returned
    T = DataFrame(t₁=categorical(["AG", "GG", "AA", "AA", "AG"]))

    t_target = TMLEEpistasis.encode(T)
    @test t_target == ["AG", "GG", "AA", "AA", "AG"]
end

@testset "Test TreatmentCountsRepr" begin
    T = DataFrame(
        RS_100=categorical(["AA", "AG", "GG", "AA", "AG"]),
        RS_10=categorical(["CC", "CG", "CG", "CC", "GG"])
        )
    
    @test TMLEEpistasis.TreatmentCountsRepr(T) ==
        "RS_100 & RS_10: AA CC 2 | AG CG 1 | AG GG 1 | GG CG 1"

end

end;

true