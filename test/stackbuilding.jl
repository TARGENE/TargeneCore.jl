module TestsStackBuilding

using Test
using BGEN
using CSV
using GenesInteraction
using MLJ
using TOML
using Distributions


GLMClassifier = @load LinearBinaryClassifier pkg=GLM verbosity=0
GLMRegressor = @load LinearRegressor pkg=GLM verbosity=0

LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity=0
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels verbosity=0
LinearRegressor = @load LinearRegressor pkg=MLJLinearModels verbosity=0
KNNRegressor = @load KNNRegressor pkg=NearestNeighborModels verbosity=0


@testset "Categorical target TMLE built from configuration file" begin
    tmle_config = joinpath("config", "tmle_categorical.toml")
    tmle =  GenesInteraction.tmle_from_toml(TOML.parsefile(tmle_config))
    @test tmle.fluctuation.glm isa GLMClassifier

    # Checking Qstack
    ## Checking Qstack.metalearner
    @test tmle.Q̅.metalearner isa LogisticClassifier
    @test tmle.Q̅.metalearner.fit_intercept == false
    ## Checking Qstack.resampling
    @test tmle.Q̅.resampling isa StratifiedCV
    @test tmle.Q̅.resampling.nfolds == 3
    ## Checking Qstack KNN models
    @test tmle.Q̅.KNNClassifier_1.K == 3
    @test tmle.Q̅.KNNClassifier_2.K == 5
    @test tmle.Q̅.KNNClassifier_3.K == 10
    ## Checking Qstack Logistic models
    @test tmle.Q̅.LogisticClassifier_1.lambda == 1.0
    @test tmle.Q̅.LogisticClassifier_1.fit_intercept == true
    @test tmle.Q̅.LogisticClassifier_2.lambda == 1.0
    @test tmle.Q̅.LogisticClassifier_2.fit_intercept == false
    @test tmle.Q̅.LogisticClassifier_3.lambda == 10
    @test tmle.Q̅.LogisticClassifier_3.fit_intercept == true
    @test tmle.Q̅.LogisticClassifier_4.lambda == 10
    @test tmle.Q̅.LogisticClassifier_4.fit_intercept == false

    # Checking Gstack
    ## Checking Gstack.metalearner
    @test tmle.G.model.metalearner isa LogisticClassifier
    @test tmle.G.model.metalearner.fit_intercept == false
    ## Checking Gstack.resampling
    @test tmle.G.model.resampling isa CV
    @test tmle.G.model.resampling.nfolds == 2
    ## Checking Gstack KNN models
    @test tmle.G.model.KNNClassifier_1.K == 3
    ## Checking Gstack Logistic models
    @test tmle.G.model.LogisticClassifier_1.lambda == 1.0
    @test tmle.G.model.LogisticClassifier_1.fit_intercept == true
    @test tmle.G.model.LogisticClassifier_2.lambda == 1.0
    @test tmle.G.model.LogisticClassifier_2.fit_intercept == false
end

@testset "Continuous target TMLE built from configuration file" begin
    tmle_config = joinpath("config", "tmle_continuous.toml")
    tmle =  GenesInteraction.tmle_from_toml(TOML.parsefile(tmle_config))
    @test tmle.fluctuation.glm isa GLMRegressor

    # Checking Qstack
    ## Checking Qstack.metalearner
    @test tmle.Q̅.metalearner isa LinearRegressor
    @test tmle.Q̅.metalearner.fit_intercept == false
    ## Checking Qstack.resampling
    @test tmle.Q̅.resampling isa CV
    @test tmle.Q̅.resampling.nfolds == 3
    ## Checking Qstack KNN models
    @test tmle.Q̅.KNNRegressor_1.K == 3
    @test tmle.Q̅.KNNRegressor_1.leafsize == 20
    @test tmle.Q̅.KNNRegressor_2.K == 5
    @test tmle.Q̅.KNNRegressor_2.leafsize == 20
    ## Checking Qstack Linear model
    @test tmle.Q̅.LinearRegressor_1.fit_intercept == true

    # Checking Gstack
    ## Checking Gstack.metalearner
    @test tmle.G.model.metalearner isa LogisticClassifier
    @test tmle.G.model.metalearner.fit_intercept == false
    ## Checking Gstack.resampling
    @test tmle.G.model.resampling isa StratifiedCV
    @test tmle.G.model.resampling.nfolds == 2
    ## Checking Gstack KNN models
    @test tmle.G.model.KNNClassifier_1.K == 3
    ## Checking Gstack Logistic models
    @test tmle.G.model.LogisticClassifier_1.lambda == 1.0
    @test tmle.G.model.LogisticClassifier_1.fit_intercept == true
    @test tmle.G.model.LogisticClassifier_2.lambda == 1.0
    @test tmle.G.model.LogisticClassifier_2.fit_intercept == false
end


end;

true

