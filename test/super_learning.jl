using GenesInteraction
using Test
using Random, Distributions
using DataFrames, CategoricalArrays
using MLJ, MLJBase


LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels
DecisionTreeClassifier = @load DecisionTreeClassifier pkg=DecisionTree


LinearRegressor = @load LinearRegressor pkg=MLJLinearModels
DecisionTreeRegressor = @load DecisionTreeRegressor pkg=DecisionTree
KNNRegressor = @load KNNRegressor pkg=NearestNeighborModels


function simulation_dataset(;n=10000, type=:classification)
    Age = rand(DiscreteUniform(30, 70), n)
    X_1 = 1 ./ (1 .+ exp.(-0.005*Age .- 0.01)) .< rand(Uniform(0,1), n)
    X_2 = 1 ./ (1 .+ exp.(-0.01*Age)) .< rand(Uniform(0,1), n)
    X = DataFrame(
        locus_1 = categorical(X_1),
        locus_2 = categorical(X_2),
        age = Age
    )

    if type == :classification
        y = 1 ./ (1 .+ exp.(-0.01*Age - 0.94*X_1 + 0.4*X_2)) .< rand(Uniform(0,1), n)
        y = categorical(y)
    else
        μ = 0.023*Age - 0.45*X_1 + 0.98*X_2
        y = μ + randn(n)*0.1
    end
    X, y
    
end


@testset "Testing the @superlearner macro" begin
    @testset "Basic usage" begin
    sl = @superlearner(LogisticClassifier(lambda=1), 
                        metalearner=LogisticClassifier(),
                        crossval=StratifiedCV(;nfolds=10, shuffle=false),
                        name=:TestSuperLearner)
    
    @test sl.logistic_classifier isa LogisticClassifier
    @test sl.metalearner isa LogisticClassifier
    @test sl.crossval isa StratifiedCV
    @test sl.crossval.nfolds == 10
    @test sl.crossval.shuffle == false
    end

    @testset "More advanced usage" begin
        sl = @superlearner(LogisticClassifier(lambda=10), 
                            LogisticClassifier(),
                            metalearner=LogisticClassifier(), 
                            name=:AdvancedSuperLearner,
                            crossval=CV(;nfolds=3, shuffle=true))
        
        @test sl.logistic_classifier isa LogisticClassifier
        @test sl.logistic_classifier.lambda == 10
        @test sl.logistic_classifier2 isa LogisticClassifier
        @test sl.metalearner isa LogisticClassifier
        @test sl.crossval isa CV
        @test sl.crossval.nfolds == 3
        @test sl.crossval.shuffle == true
        @test library(sl) == (:logistic_classifier, :logistic_classifier2)
        end

end


@testset "Testing fit! method of the superlearner" begin
    @testset "With logistic regression simulation dataset" begin
        Random.seed!(123)
        X, y = simulation_dataset(;n=10000)
        pipe = @pipeline OneHotEncoder(;drop_last=true) @superlearner(LogisticClassifier(lambda=1),
                        ConstantClassifier(),
                        KNNClassifier(),
                        metalearner=LogisticClassifier(),
                        name=:Sim1SuperLearner,
                        crossval=StratifiedCV(;nfolds=3, shuffle=false)
                ) name=MyPipeline
        mach = machine(pipe, X, y)
        MLJ.fit!(mach)
        fp = fitted_params(mach)

        # Retrieve metalearner coefficients and check that the true model
        # receives the highest coefficient
        metalearner_coefs = fp.sim1_super_learner.metalearner[end].coefs
        @test metalearner_coefs[1].second > 3
        for i in 2:3
            @test metalearner_coefs[i].second < 2
        end

        # Check coefficients of the true model, the true model definitions generates
        # data according to X1,X2=true and X2=true while the model fits X1,X2=false 
        lr_coefs = fp.sim1_super_learner.logistic_classifier[end-1].coefs
        for pair in lr_coefs
            if pair.first == :locus_1__false
                @test pair.second ≈ 0.94 atol=0.006
            elseif pair.first == :locus_2__false
                @test pair.second ≈ -0.4 atol=0.006
            elseif pair.first == :age
                @test pair.second ≈ -0.01 atol=0.006
            end
        end
    end

    @testset "With linear regression simulation dataset" begin
        Random.seed!(123)
        X, y = simulation_dataset(;n=10000, type=:regression)
        pipe = @pipeline OneHotEncoder(;drop_last=true) @superlearner(LinearRegressor(),
                        DecisionTreeRegressor(),
                        KNNRegressor(),
                        metalearner=LinearRegressor(),
                        name=:SimRegSuperLearner,
                        crossval=CV(;nfolds=3, shuffle=false)
                ) name=MyRegPipeline
        mach = machine(pipe, X, y)
        MLJ.fit!(mach)

        fp = fitted_params(mach)

        # Retrieve metalearner coefficients and check that the true model
        # receives the highest coefficient
        metalearner_coefs = fp.sim_reg_super_learner.metalearner[end].coefs
        @test metalearner_coefs[1].second > 1
        for i in 2:3
            @test abs(metalearner_coefs[i].second) < 0.15
        end
        
        # Check coefficients of the true model, the true model definitions generates
        # data according to X1,X2=true and X2=true while the model fits X1,X2=false 
        lr_coefs = fp.sim_reg_super_learner.linear_regressor[end-1].coefs
        for pair in lr_coefs
            if pair.first == :locus_1__false
                @test pair.second ≈ 0.45 atol=0.006
            elseif pair.first == :locus_2__false
                @test pair.second ≈ -0.98 atol=0.006
            elseif pair.first == :age
                @test pair.second ≈ 0.023 atol=0.006
            end
        end
    end
end

@testset "Testing helper functions" begin
    y = source()

    # Checking the return type is a Node
    # Maybe test further the associated operation sometime? Need to know about theinternals.
    sl_cv = @superlearner(LogisticClassifier(), 
                        metalearner=LogisticClassifier(), 
                        name=:SuperLearnerCV,
                        crossval=CV())

    @test GenesInteraction.getfolds(y, sl_cv, 1000) isa Source

    sl_scv = @superlearner(LogisticClassifier(), 
                        metalearner=LogisticClassifier(), 
                        name=:SuperLearnerSCV,
                        crossval=StratifiedCV())
    @test GenesInteraction.getfolds(y, sl_scv, 1000) isa Node

    @test library(sl_scv) == (:logistic_classifier,)
end