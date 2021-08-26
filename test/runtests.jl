using Test
using GenesInteraction
using MLJ
using Distributions
using TOML
using CSV
using DataFrames
using CategoricalArrays

genotypefile = joinpath("data", "mouse")
phenotypefile = joinpath("data", "pheno_10_causals.txt")
confoundersfile = joinpath("data", "confouders.csv")
snpfile = joinpath("data", "query_snps.csv")

LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity = 0
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels verbosity = 0
LinearRegressor = @load LinearRegressor pkg=MLJLinearModels verbosity = 0
KNNRegressor = @load KNNRegressor pkg=NearestNeighborModels verbosity = 0


@testset "Categorical target TMLE built from configuration file" begin
    tmle_config = joinpath("config", "categorical.toml")
    tmle =  GenesInteraction.tmle_from_toml(TOML.parsefile(tmle_config))
    @test tmle.fluctuation_family isa Bernoulli

    # Checking Qstack
    Qstack = tmle.target_cond_expectation_estimator
    ## Checking Qstack.metalearner
    @test Qstack.metalearner isa LogisticClassifier
    @test Qstack.metalearner.fit_intercept == false
    ## Checking Qstack.resampling
    @test Qstack.resampling isa StratifiedCV
    @test Qstack.resampling.nfolds == 3
    ## Checking Qstack KNN models
    @test Qstack.KNNClassifier_1.K == 3
    @test Qstack.KNNClassifier_2.K == 5
    @test Qstack.KNNClassifier_3.K == 10
    ## Checking Qstack Logistic models
    @test Qstack.LogisticClassifier_1.lambda == 1.0
    @test Qstack.LogisticClassifier_1.fit_intercept == true
    @test Qstack.LogisticClassifier_2.lambda == 1.0
    @test Qstack.LogisticClassifier_2.fit_intercept == false
    @test Qstack.LogisticClassifier_3.lambda == 10
    @test Qstack.LogisticClassifier_3.fit_intercept == true
    @test Qstack.LogisticClassifier_4.lambda == 10
    @test Qstack.LogisticClassifier_4.fit_intercept == false

    # Checking Gstack
    Gstack = tmle.treatment_cond_likelihood_estimator
    ## Checking Gstack.metalearner
    @test Gstack.metalearner isa LogisticClassifier
    @test Gstack.metalearner.fit_intercept == false
    ## Checking Gstack.resampling
    @test Gstack.resampling isa CV
    @test Gstack.resampling.nfolds == 2
    ## Checking Gstack KNN models
    @test Gstack.KNNClassifier_1.K == 3
    ## Checking Gstack Logistic models
    @test Gstack.LogisticClassifier_1.lambda == 1.0
    @test Gstack.LogisticClassifier_1.fit_intercept == true
    @test Gstack.LogisticClassifier_2.lambda == 1.0
    @test Gstack.LogisticClassifier_2.fit_intercept == false
end

@testset "Continuous target TMLE built from configuration file" begin
    tmle_config = joinpath("config", "continuous.toml")
    tmle =  GenesInteraction.tmle_from_toml(TOML.parsefile(tmle_config))
    @test tmle.fluctuation_family isa Normal

    # Checking Qstack
    Qstack = tmle.target_cond_expectation_estimator
    ## Checking Qstack.metalearner
    @test Qstack.metalearner isa LinearRegressor
    @test Qstack.metalearner.fit_intercept == false
    ## Checking Qstack.resampling
    @test Qstack.resampling isa CV
    @test Qstack.resampling.nfolds == 3
    ## Checking Qstack KNN models
    @test Qstack.KNNRegressor_1.K == 3
    @test Qstack.KNNRegressor_1.leafsize == 20
    @test Qstack.KNNRegressor_2.K == 5
    @test Qstack.KNNRegressor_2.leafsize == 20
    ## Checking Qstack Linear model
    @test Qstack.LinearRegressor_1.fit_intercept == true

    # Checking Gstack
    Gstack = tmle.treatment_cond_likelihood_estimator
    ## Checking Gstack.metalearner
    @test Gstack.metalearner isa LogisticClassifier
    @test Gstack.metalearner.fit_intercept == false
    ## Checking Gstack.resampling
    @test Gstack.resampling isa StratifiedCV
    @test Gstack.resampling.nfolds == 2
    ## Checking Gstack KNN models
    @test Gstack.KNNClassifier_1.K == 3
    ## Checking Gstack Logistic models
    @test Gstack.LogisticClassifier_1.lambda == 1.0
    @test Gstack.LogisticClassifier_1.fit_intercept == true
    @test Gstack.LogisticClassifier_2.lambda == 1.0
    @test Gstack.LogisticClassifier_2.fit_intercept == false
end


@testset "prepareTreatment" begin
    snpqueries = CSV.File(snpfile) |> DataFrame
    T = [2 0;
        0 0;
        2 3;
        0 3]
    T = GenesInteraction.prepareTreatment(T, first(snpqueries))
    @test Tables.getcolumn(T, :x1) == [1, 0, 1, 0]
    @test Tables.getcolumn(T, :x2) == [0, 0, 1, 1]
end

@testset "prepareTarget" begin
    # categorical config
    y = [true, false, true]
    config = TOML.parsefile(joinpath("config", "categorical.toml"))
    y = GenesInteraction.prepareTarget(y, config)
    @test y == CategoricalArray{Bool,1,UInt32}([true, false, true])
    # continuous config
    y = [1, 0, 1]
    config = TOML.parsefile(joinpath("config", "continuous.toml"))
    @test GenesInteraction.prepareTarget(y, config) == y
end

@testset "Test filtering_idx" begin
    snparray = [3 2 3 0;
                0 2 0 1;
                1 1 1 1;
                3 0 0 1;
                2 3 3 0;
                0 0 2 2]
    snpqueries = CSV.File(snpfile) |> DataFrame
    snpquery = first(snpqueries)
    snp_idx = Dict("rs3683945" => 2, "rs3706761"=> 1)
    indices, columns = GenesInteraction.filtering_idx(snparray, snpquery, snp_idx)
    @test columns == [2, 1]
    @test indices == [1, 2, 4, 6]
end

@testset "Epistasis estimation run" begin

    tmle_config = joinpath("config", "continuous.toml")
    outfile = "tmleepistasis_results.csv"
    # results = @enter tmleepistasis(genotypefile, 
    #                 phenotypefile, 
    #                 confoundersfile, 
    #                 snpfile,
    #                 tmle_config,
    #                 outfile)

end