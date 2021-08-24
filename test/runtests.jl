using Test
using GenesInteraction


genfile = joinpath("data", "mouse")
phenotypefile = joinpath("data", "pheno_10_causals.txt")
confoundersfile = joinpath("data", "confouders.csv")
snpfile = joinpath("data", "query_snps.csv")
tmle_config = joinpath("config", "categorical_test.toml")

@testset "Testing TMLE from configuratin file" begin
    @enter GenesInteraction.tmle_from_toml(tmle_config)
end

results = runepistatis(genfile, phenotypefile, confoundersfile, snpfile)
