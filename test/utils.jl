module TestUtils

using Test
using GenesInteraction
using MLJ
using TOML

@testset "Test retrieve coefficients" begin
    estimatorfile = joinpath("config", "tmle_config.toml")
    tmle_config = TOML.parsefile(estimatorfile)
    queries = ["query_1" => (t₁=[1, 0], t₂=[1, 0]),]
    tmle = GenesInteraction.tmles_from_toml(tmle_config, queries)

    n = 1000
    T = (t₁=categorical(rand([1, 0], n)), t₂=categorical(rand([1, 0], n)))
    W = (w₁=rand(n),)
    y = rand(n)
    mach = machine(tmle["continuous"], T, W, y)
    fit!(mach, verbosity=0)

    coefs = GenesInteraction.Qstack_coefs(mach)
    @test coefs isa Vector
    @test length(coefs) == 4

    repr = GenesInteraction.repr_Qstack_coefs(mach)
    @test repr isa String
    str_reprs = split(repr, ", ")
    modelstr = join([split(v, " => ")[1] for v in str_reprs], ",")
    @test modelstr == "HALRegressor,XGBoostRegressor,XGBoostRegressor,InteractionLMRegressor"

end

@testset "Test parse_queries" begin
    build_query_file()
    queries = GenesInteraction.parse_queries(queryfile)
    @test queries == [
        "QUERY_1" => (RSID_10 = ["AG", "GG"], RSID_100 = ["AG", "GG"]),
        "QUERY_2" => (RSID_10 = ["AG", "GG"], RSID_100 = ["AA", "GG"])
    ]
    rm(queryfile)

    @test GenesInteraction.querystring(queries[1][2]) == "RSID_10: AG -> GG & RSID_100: AG -> GG"
end

@testset "Test phenotypes parsing" begin
    allnames = GenesInteraction.phenotypesnames(phenotypefile)
    @test allnames == [:categorical_phenotype, :continuous_phenotype]

    @test GenesInteraction.phenotypes_list(nothing, [:categorical_phenotype], allnames) == [:continuous_phenotype]
    @test GenesInteraction.phenotypes_list(nothing, [], allnames) == allnames

    @test GenesInteraction.phenotypes_list("data/phen_list_1.csv", [], allnames) == [:continuous_phenotype]
    @test GenesInteraction.phenotypes_list("data/phen_list_2.csv", [:continuous_phenotype], allnames) == [:categorical_phenotype]

    outfile = "temp_results.csv"
    @test GenesInteraction.init_or_retrieve_results(outfile) == Set()
    CSV.write(outfile, DataFrame(PHENOTYPE=["toto", "tata", "toto"]))
    @test GenesInteraction.init_or_retrieve_results(outfile) == Set([:toto, :tata])
    rm(outfile)

end

end;

true