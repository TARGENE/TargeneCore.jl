using GenesInteraction
using Test
using Random, Distributions
using DataFrames, CategoricalArrays
using MLJ

import MLJLinearModels.LogisticClassifier
import DecisionTree.DecisionTreeClassifier
import EvoTrees.EvoTreeClassifier
import NearestNeighborModels.KNNClassifier


function simulation_dataset(;n=10000)
    Age = rand(DiscreteUniform(30, 70), n)
    X_1 = 1 ./ (1 .+ exp.(-0.005*Age .- 0.01)) .< rand(Uniform(0,1), n)
    X_2 = 1 ./ (1 .+ exp.(-0.01*Age)) .< rand(Uniform(0,1), n)
    X = DataFrame(
        locus_1 = categorical(X_1),
        locus_2 = categorical(X_2),
        age = Age
    )

    y = 1 ./ (1 .+ exp.(-0.01*Age - 0.94*X_1 + 0.4*X_2)) .< rand(Uniform(0,1), n)
    X, categorical(y)
end


@testset "Testing the @superlearner macro" begin
    @testset "Basic usage" begin
    sl = @superlearner(LogisticClassifier(lambda=1), 
                        metalearner=LogisticClassifier(), 
                        name=:TestSuperLearner)
    
    @test sl.logistic_classifier isa LogisticClassifier
    @test sl.metalearner isa LogisticClassifier
    @test sl.nfolds == 10
    @test sl.shuffle == false
    end

    @testset "More advanced usage" begin
        sl = @superlearner(LogisticClassifier(lambda=10), 
                            LogisticClassifier(),
                            metalearner=LogisticClassifier(), 
                            name=:AdvancedSuperLearner,
                            nfolds=3,
                            shuffle=true)
        
        @test sl.logistic_classifier isa LogisticClassifier
        @test sl.logistic_classifier.lambda == 10
        @test sl.logistic_classifier2 isa LogisticClassifier
        @test sl.metalearner isa LogisticClassifier
        @test sl.nfolds == 3
        @test sl.shuffle == true
        end

end


@testset "Testing fit! method of the superlearner" begin
    @testset "With logistic regression simulation dataset" begin
        Random.seed!(123)
        X, y = simulation_dataset(;n=10000)
        pipe = @pipeline OneHotEncoder(;drop_last=true) @superlearner(LogisticClassifier(lambda=1),
                        ConstantClassifier(),
                        KNNClassifier(),
                        EvoTreeClassifier(),
                        metalearner=LogisticClassifier(),
                        name=:Sim1SuperLearner,
                        nfolds=10,
                        shuffle=false
                ) name=MyPipeline
        mach = machine(pipe, X, y)
        MLJ.fit!(mach)
        fp = fitted_params(mach)

        # Retrieve metalearner coefficients and check that the true model
        # receives the highest coefficient
        metalearner_coefs = fp.sim1_super_learner.metalearner[end].coefs
        @test metalearner_coefs[1].second > 3
        for i in 2:4
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

    @testset "Emphasizing the problem if the same model is used multiple times"  begin
        Random.seed!(123)
        X, y = simulation_dataset(;n=10000)
        pipe = @pipeline OneHotEncoder(;drop_last=true) @superlearner(LogisticClassifier(lambda=1),
                        ConstantClassifier(),
                        KNNClassifier(),
                        EvoTreeClassifier(),
                        metalearner=LogisticClassifier(),
                        name=:PbSuperLearner,
                        nfolds=10,
                        shuffle=false
                ) name=MyPbPipeline
        mach = machine(pipe, X, y)
        MLJ.fit!(mach)
        fp = fitted_params(mach)

        # The metalearner has only been fitted once and the logistic regression 11 times
        # I think because they share the same model, the current fitted_params method can't 
        # make a difference between the 2.
        # What is the unit of identification?
        @test fp.pb_super_learner.metalearner == fp.pb_super_learner.logistic_classifier
    end
end