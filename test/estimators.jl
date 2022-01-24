module TestsStackBuilding

using Test
using BGEN
using CSV
using TMLEEpistasis
using MLJ
using TOML
using TMLE
using Distributions

include("helper_fns.jl")

GLMClassifier = @load LinearBinaryClassifier pkg=GLM verbosity=0
LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity=0

GLMRegressor = @load LinearRegressor pkg=GLM verbosity=0
LinearRegressor = @load LinearRegressor pkg=MLJLinearModels verbosity=0


@testset "Categorical target TMLE built from configuration file" begin
    tmle_config = joinpath("config", "tmle_config.toml")
    build_query_file()
    queries = TMLEEpistasis.parse_queries(queryfile)
    tmles =  TMLEEpistasis.estimators_from_toml(TOML.parsefile(tmle_config), queries)
    # Test binary target TMLE's Qstack
    tmle = tmles["binary"]
    @test tmle.F isa GLMClassifier
    ## Checking Qstack.metalearner
    @test tmle.Q̅.metalearner isa LogisticClassifier
    @test tmle.Q̅.metalearner.fit_intercept == false
    ## Checking Qstack.resampling
    @test tmle.Q̅.resampling isa StratifiedCV
    @test tmle.Q̅.resampling.nfolds == 2
    ## Checking Qstack XGBoost models
    @test tmle.Q̅.XGBoostClassifier_1.num_round == 10
    ## Checking Qstack  Interaction Logistic models
    @test tmle.Q̅.InteractionLMClassifier_1 isa TMLEEpistasis.InteractionLMClassifier
    @test tmle.Q̅.InteractionLMClassifier_1.interaction_transformer.column_pattern == r"^RS_"
    ## Checking Qstack HAL model
    @test tmle.Q̅.HALClassifier_1.lambda == 10
    @test tmle.Q̅.HALClassifier_1.smoothness_orders == 1
    @test tmle.Q̅.HALClassifier_1.cv_select == false
    @test tmle.Q̅.HALClassifier_1.num_knots == [10, 5]

    # Test binary target TMLE Qstack
    tmle = tmles["continuous"]
    @test tmle.F isa GLMRegressor
    ## Checking Qstack.metalearner
    @test tmle.Q̅.metalearner isa LinearRegressor
    @test tmle.Q̅.metalearner.fit_intercept == false

    ## Checking Qstack.resampling
    @test tmle.Q̅.resampling isa CV
    @test tmle.Q̅.resampling.nfolds == 2
    ## Checking Qstack XGBoost models
    @test tmle.Q̅.XGBoostRegressor_1.num_round == 10
    @test tmle.Q̅.XGBoostRegressor_2.num_round == 20
    ## Checking Qstack Interaction Linear model
    @test tmle.Q̅.InteractionLMRegressor_1.interaction_transformer.column_pattern == r"^RS_"
    ## Checking Qstack HAL model
    @test tmle.Q̅.HALRegressor_1.lambda == 10
    @test tmle.Q̅.HALRegressor_1.smoothness_orders == 1
    @test tmle.Q̅.HALRegressor_1.cv_select == false
    @test tmle.Q̅.HALRegressor_1.num_knots == [10, 5]
    
    # Both TMLE have the same G Stack
    for (type, tmle) in tmles
        @test tmle.queries == (
            (RSID_10 = ["AG", "GG"], RSID_100 = ["AG", "GG"]), 
            (RSID_10 = ["AG", "GG"], RSID_100 = ["AA", "GG"])
            )
        # Checking Gstack
        @test tmle.G isa FullCategoricalJoint 
        ## Checking Gstack.metalearner
        @test tmle.G.model.metalearner isa LogisticClassifier
        @test tmle.G.model.metalearner.fit_intercept == false
        ## Checking Gstack.resampling
        @test tmle.G.model.resampling isa StratifiedCV
        @test tmle.G.model.resampling.nfolds == 2
        ## Checking Gstack models
        @test tmle.G.model.LogisticClassifier_1.lambda == 1.0
        @test tmle.G.model.XGBoostClassifier_1.num_round == 10
    end
    rm(queryfile)
end

@testset "Test standard ATE TMLE build" begin
    tmle_config = joinpath("config", "tmle_config.toml")
    build_ate_query_file()
    queries = TMLEEpistasis.parse_queries(queryfile)
    tmles =  TMLEEpistasis.estimators_from_toml(TOML.parsefile(tmle_config), queries)
    for (type, tmle) in tmles
        @test tmle.queries == (
            (RSID_10 = ["AG", "GG"],), 
            (RSID_10 = ["AA", "GG"],)
            )
        # Checking Gstack
        @test tmle.G isa Stack
    end
    rm(queryfile)
end



end;

true

