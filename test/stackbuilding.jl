module TestsStackBuilding

using Test
using BGEN
using CSV
using GenesInteraction
using MLJ
using TOML
using Distributions


GLMClassifier = @load LinearBinaryClassifier pkg=GLM verbosity=0
SKLogisticClassifier = @load LogisticClassifier pkg=ScikitLearn verbosity=0
LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity=0
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels verbosity=0

GLMRegressor = @load LinearRegressor pkg=GLM verbosity=0
SKLinearRegressor = @load LinearRegressor pkg=ScikitLearn verbosity=0
LinearRegressor = @load LinearRegressor pkg=MLJLinearModels verbosity=0
KNNRegressor = @load KNNRegressor pkg=NearestNeighborModels verbosity=0


@testset "Categorical target TMLE built from configuration file" begin
    tmle_config = joinpath("config", "tmle_config.toml")

    tmles =  GenesInteraction.tmles_from_toml(TOML.parsefile(tmle_config))
    # Test binary target TMLE's Qstack
    tmle = tmles["binary"]
    @test tmle.F.glm isa GLMClassifier
    ## Checking Qstack.metalearner
    @test tmle.Q̅.metalearner isa SKLogisticClassifier
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
    ## Checking Qstack HAL model
    @test tmle.Q̅.HALClassifier_1.lambda == 10
    @test tmle.Q̅.HALClassifier_1.smoothness_orders == 2
    @test tmle.Q̅.HALClassifier_1.cv_select == false
    @test tmle.Q̅.HALClassifier_1.num_knots == [10, 5]

    # Test binary target TMLE Qstack
    tmle = tmles["continuous"]
    @test tmle.F.glm isa GLMRegressor
    ## Checking Qstack.metalearner
    @test tmle.Q̅.metalearner isa SKLinearRegressor
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
    ## Checking Qstack HAL model
    @test tmle.Q̅.HALRegressor_1.lambda == 10
    @test tmle.Q̅.HALRegressor_1.smoothness_orders == 2
    @test tmle.Q̅.HALRegressor_1.cv_select == false
    @test tmle.Q̅.HALRegressor_1.num_knots == [10, 5]
    
    # Both TMLE have the same G Stack
    for (type, tmle) in tmles
        # Checking Gstack
        ## Checking Gstack.metalearner
        @test tmle.G.model.metalearner isa SKLogisticClassifier
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
end


end;

true

