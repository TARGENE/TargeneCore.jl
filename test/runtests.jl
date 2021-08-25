using Test
using GenesInteraction
using MLJ
using Distributions


LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity = 0
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels verbosity = 0

genfile = joinpath("data", "mouse")
phenotypefile = joinpath("data", "pheno_10_causals.txt")
confoundersfile = joinpath("data", "confouders.csv")
snpfile = joinpath("data", "query_snps.csv")
tmle_config = joinpath("config", "categorical_test.toml")


@testset "Testing TMLE from configuratin file" begin
    tmle =  GenesInteraction.tmle_from_toml(tmle_config)
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
    ## Checking Qstack.metalearner
    @test Gstack.metalearner isa LogisticClassifier
    @test Gstack.metalearner.fit_intercept == false
    ## Checking Qstack.resampling
    @test Gstack.resampling isa CV
    @test Gstack.resampling.nfolds == 2
    ## Checking Qstack KNN models
    @test Gstack.KNNClassifier_1.K == 3
    ## Checking Qstack Logistic models
    @test Gstack.LogisticClassifier_1.lambda == 1.0
    @test Gstack.LogisticClassifier_1.fit_intercept == true
    @test Gstack.LogisticClassifier_2.lambda == 1.0
    @test Gstack.LogisticClassifier_2.fit_intercept == false
end

#results = runepistatis(genfile, phenotypefile, confoundersfile, snpfile)
