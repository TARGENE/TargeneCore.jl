module TestFromActors

using Test
using CSV
using DataFrames
using TargeneCore
using YAML
using TMLE

include("test_utils.jl")


#####################################################################
###############               UNIT TESTS              ###############
#####################################################################

@testset "Test combine_by_bqtl" begin
    bQTLS = DataFrame(
        ID = ["rs1", "rs2", "rs3"],
        CHR = [1, 2, 3],
        DESCRIPTION = ["VDR", "VDR", "VDR"]
    )
    trans_actors = [
        DataFrame(
            ID = ["rs4", "rs5"],
            CHR = [1, 12],
            DESCRIPTION = ["VitD-QTL", "VitD-QTL"]
    ),
        DataFrame(
            ID = ["rs6", "rs7", "rs8"],
            CHR = [12, 9, 8],
            DESCRIPTION = ["RXR-QTL", "RXR-QTL", "RXR-QTL"]
        )
    ]
    envs = DataFrame(ID=["sex"], DESCRIPTION=["Env"])
    # order 1
    combinations = TargeneCore.combine_by_bqtl(bQTLS, trans_actors, nothing, 1)
    @test combinations == [[:rs1], [:rs2], [:rs3]]
    # order 2, no environmental variable
    combinations = TargeneCore.combine_by_bqtl(bQTLS, trans_actors, nothing, 2)
    @test combinations == [[:rs1, :rs4],
                           [:rs1, :rs5],
                           [:rs2, :rs4],
                           [:rs2, :rs5],
                           [:rs3, :rs4],
                           [:rs3, :rs5],
                           [:rs1, :rs6],
                           [:rs1, :rs7],
                           [:rs1, :rs8],
                           [:rs2, :rs6],
                           [:rs2, :rs7],
                           [:rs2, :rs8],
                           [:rs3, :rs6],
                           [:rs3, :rs7],
                           [:rs3, :rs8]]
    # order 2, only one trans-actors dataset
    combinations = TargeneCore.combine_by_bqtl(bQTLS, [trans_actors[1]], nothing, 2)
    @test combinations == [[:rs1, :rs4],
                           [:rs1, :rs5],
                           [:rs2, :rs4],
                           [:rs2, :rs5],
                           [:rs3, :rs4],
                           [:rs3, :rs5]]
    # order 3, no environmental variable
    combinations = TargeneCore.combine_by_bqtl(bQTLS, trans_actors, nothing, 3)
    @test combinations == [[:rs1, :rs4, :rs6],
                           [:rs1, :rs4, :rs7],
                           [:rs1, :rs4, :rs8],
                           [:rs1, :rs5, :rs6],
                           [:rs1, :rs5, :rs7],
                           [:rs1, :rs5, :rs8],
                           [:rs2, :rs4, :rs6],
                           [:rs2, :rs4, :rs7],
                           [:rs2, :rs4, :rs8],
                           [:rs2, :rs5, :rs6],
                           [:rs2, :rs5, :rs7],
                           [:rs2, :rs5, :rs8],
                           [:rs3, :rs4, :rs6],
                           [:rs3, :rs4, :rs7],
                           [:rs3, :rs4, :rs8],
                           [:rs3, :rs5, :rs6],
                           [:rs3, :rs5, :rs7],
                           [:rs3, :rs5, :rs8]]
    
    # order 2, environmental variable
    combinations = TargeneCore.combine_by_bqtl(bQTLS, trans_actors, envs, 2)
    @test combinations == [[:rs1, :rs4],
                           [:rs1, :rs5],
                           [:rs2, :rs4],
                           [:rs2, :rs5],
                           [:rs3, :rs4],
                           [:rs3, :rs5],
                           [:rs1, :rs6],
                           [:rs1, :rs7],
                           [:rs1, :rs8],
                           [:rs2, :rs6],
                           [:rs2, :rs7],
                           [:rs2, :rs8],
                           [:rs3, :rs6],
                           [:rs3, :rs7],
                           [:rs3, :rs8],
                           [:rs1, :sex],
                           [:rs2, :sex],
                           [:rs3, :sex]]
    # order 3, only one trans-actors dataset
    combinations = TargeneCore.combine_by_bqtl(bQTLS, [trans_actors[1]], envs, 3)
    @test combinations == [[:rs1, :rs4, :sex],
                           [:rs1, :rs5, :sex],
                           [:rs2, :rs4, :sex],
                           [:rs2, :rs5, :sex],
                           [:rs3, :rs4, :sex],
                           [:rs3, :rs5, :sex]]
    # order 2, no trans-actors dataset
    combinations = TargeneCore.combine_by_bqtl(bQTLS, nothing, envs, 2)
    @test combinations == [[:rs1, :sex],
                           [:rs2, :sex],
                           [:rs3, :sex]]
end

@testset "Test control_case_settings" begin
    data = DataFrame(
        T1 = [0, 1, 2],
        T2 = ["AC", "CC", "AA"],
        T3 = [0, 1, 1] 
    )
    # One Treatment ATE
    treatments = [:T1]
    settings = collect(TargeneCore.control_case_settings(ATE, treatments, data))
    @test settings == [
        (T1 = (case = 0, control = 1),),
        (T1 = (case = 0, control = 2),),
        (T1 = (case = 1, control = 2),)
    ]
    # Two treatments IATEs
    treatments = [:T1, :T2]
    settings = collect(TargeneCore.control_case_settings(IATE, treatments, data))
    expected_settings = [
        (T1 = (case = 0, control = 1), T2 = (case = "AA", control = "AC")),
        (T1 = (case = 0, control = 2), T2 = (case = "AA", control = "AC")),
        (T1 = (case = 1, control = 2), T2 = (case = "AA", control = "AC")),
        (T1 = (case = 0, control = 1), T2 = (case = "AA", control = "CC")),
        (T1 = (case = 0, control = 2), T2 = (case = "AA", control = "CC")),
        (T1 = (case = 1, control = 2), T2 = (case = "AA", control = "CC")),
        (T1 = (case = 0, control = 1), T2 = (case = "AC", control = "CC")),
        (T1 = (case = 0, control = 2), T2 = (case = "AC", control = "CC")),
        (T1 = (case = 1, control = 2), T2 = (case = "AC", control = "CC"))
    ]
    @test length(settings) == length(expected_settings)
    for index in eachindex(expected_settings)
        @test expected_settings[index] == settings[index]
    end

    # Two treatments ATEs
    settings = collect(TargeneCore.control_case_settings(ATE, treatments, data))
    expected_settings = [
        (T1 = (case = 0, control = 1), T2 = (case = "AA", control = "AA")),
        (T1 = (case = 0, control = 2), T2 = (case = "AA", control = "AA")),
        (T1 = (case = 1, control = 2), T2 = (case = "AA", control = "AA")),
        (T1 = (case = 0, control = 1), T2 = (case = "AC", control = "AC")),
        (T1 = (case = 0, control = 2), T2 = (case = "AC", control = "AC")),
        (T1 = (case = 1, control = 2), T2 = (case = "AC", control = "AC")),
        (T1 = (case = 0, control = 1), T2 = (case = "CC", control = "CC")),
        (T1 = (case = 0, control = 2), T2 = (case = "CC", control = "CC")),
        (T1 = (case = 1, control = 2), T2 = (case = "CC", control = "CC"))
    ]
    @test length(settings) == length(expected_settings)
    for index in eachindex(expected_settings)
        @test expected_settings[index] == settings[index]
    end

    # Three treatments IATEs
    treatments = [:T1, :T2, :T3]
    settings = collect(TargeneCore.control_case_settings(IATE, treatments, data))
    expected_settings = [
        (T1 = (case = 0, control = 1), T2 = (case = "AA", control = "AC"), T3 = (case = 0, control = 1)),
        (T1 = (case = 0, control = 2), T2 = (case = "AA", control = "AC"), T3 = (case = 0, control = 1)),
        (T1 = (case = 1, control = 2), T2 = (case = "AA", control = "AC"), T3 = (case = 0, control = 1)),
        (T1 = (case = 0, control = 1), T2 = (case = "AA", control = "CC"), T3 = (case = 0, control = 1)),
        (T1 = (case = 0, control = 2), T2 = (case = "AA", control = "CC"), T3 = (case = 0, control = 1)),
        (T1 = (case = 1, control = 2), T2 = (case = "AA", control = "CC"), T3 = (case = 0, control = 1)),
        (T1 = (case = 0, control = 1), T2 = (case = "AC", control = "CC"), T3 = (case = 0, control = 1)),
        (T1 = (case = 0, control = 2), T2 = (case = "AC", control = "CC"), T3 = (case = 0, control = 1)),
        (T1 = (case = 1, control = 2), T2 = (case = "AC", control = "CC"), T3 = (case = 0, control = 1))
    ]
    @test length(settings) == length(expected_settings)
    for index in eachindex(expected_settings)
        @test expected_settings[index] == settings[index]
    end

    # Three treatments ATEs
    settings = collect(TargeneCore.control_case_settings(ATE, treatments, data))
    expected_settings = [
        (T1 = (case = 0, control = 1), T2 = (case = "AA", control = "AA"), T3 = (case = 0, control = 0)),
        (T1 = (case = 0, control = 2), T2 = (case = "AA", control = "AA"), T3 = (case = 0, control = 0)),
        (T1 = (case = 1, control = 2), T2 = (case = "AA", control = "AA"), T3 = (case = 0, control = 0)),
        (T1 = (case = 0, control = 1), T2 = (case = "AC", control = "AC"), T3 = (case = 0, control = 0)),
        (T1 = (case = 0, control = 2), T2 = (case = "AC", control = "AC"), T3 = (case = 0, control = 0)),
        (T1 = (case = 1, control = 2), T2 = (case = "AC", control = "AC"), T3 = (case = 0, control = 0)),
        (T1 = (case = 0, control = 1), T2 = (case = "CC", control = "CC"), T3 = (case = 0, control = 0)),
        (T1 = (case = 0, control = 2), T2 = (case = "CC", control = "CC"), T3 = (case = 0, control = 0)),
        (T1 = (case = 1, control = 2), T2 = (case = "CC", control = "CC"), T3 = (case = 0, control = 0)),
        (T1 = (case = 0, control = 1), T2 = (case = "AA", control = "AA"), T3 = (case = 1, control = 1)),
        (T1 = (case = 0, control = 2), T2 = (case = "AA", control = "AA"), T3 = (case = 1, control = 1)),
        (T1 = (case = 1, control = 2), T2 = (case = "AA", control = "AA"), T3 = (case = 1, control = 1)),
        (T1 = (case = 0, control = 1), T2 = (case = "AC", control = "AC"), T3 = (case = 1, control = 1)),
        (T1 = (case = 0, control = 2), T2 = (case = "AC", control = "AC"), T3 = (case = 1, control = 1)),
        (T1 = (case = 1, control = 2), T2 = (case = "AC", control = "AC"), T3 = (case = 1, control = 1)),
        (T1 = (case = 0, control = 1), T2 = (case = "CC", control = "CC"), T3 = (case = 1, control = 1)),
        (T1 = (case = 0, control = 2), T2 = (case = "CC", control = "CC"), T3 = (case = 1, control = 1)),
        (T1 = (case = 1, control = 2), T2 = (case = "CC", control = "CC"), T3 = (case = 1, control = 1))
    ]
    @test length(settings) == length(expected_settings)
    for index in eachindex(settings)
        @test expected_settings[index] == settings[index]
    end
end

@testset "Test treatments_from_actors" begin
    # At least two type of actors should be specified
    @test_throws ArgumentError TargeneCore.treatments_from_actors(nothing, nothing, 1)
    @test_throws ArgumentError TargeneCore.treatments_from_actors(1, nothing, nothing)
    @test_throws ArgumentError TargeneCore.treatments_from_actors(nothing, 1, nothing)

    bqtl_file = joinpath("data", "bqtls.csv")
    trans_actors_prefix = joinpath("data", "trans_actors_1.csv")
    env_file = joinpath("data", "extra_treatments.txt")
    # bqtls and trans_actors
    bqtls, transactors, envs = TargeneCore.treatments_from_actors(bqtl_file, nothing, trans_actors_prefix)
    @test bqtls isa DataFrame
    @test envs isa Nothing
    @test transactors isa Vector{DataFrame}
    @test size(transactors, 1) == 1

    # bqtls and env
    bqtls, transactors, envs = TargeneCore.treatments_from_actors(bqtl_file, env_file, nothing)
    @test bqtls isa DataFrame
    @test envs == ["TREAT_1"]
    @test transactors isa Nothing

    # trans actors and env
    bqtls, transactors, envs = TargeneCore.treatments_from_actors(nothing, env_file, trans_actors_prefix)
    @test bqtls isa Nothing
    @test envs == ["TREAT_1"]
    @test transactors isa Vector{DataFrame}
    @test size(transactors, 1) == 1
end

#####################################################################
###############           END-TO-END TESTS            ###############
#####################################################################


@testset "Test tmle_inputs from-actors: scenario 1" begin
    # Scenario:
    # - Trans-actors
    # - Extra Treatment
    # - Extra Covariates
    # - Order 1,2
    parsed_args = Dict(
        "from-actors" => Dict{String, Any}(
            "bqtls" => joinpath("data", "bqtls.csv"), 
            "trans-actors-prefix" => joinpath("data", "trans_actors_1.csv"),
            "extra-covariates" => joinpath("data", "extra_covariates.txt"),
            "extra-treatments" => joinpath("data", "extra_treatments.txt"),
            "extra-confounders" => nothing,
            "orders" => "1,2",
            "genotypes-as-int" => true
            ),
        "traits" => joinpath("data", "traits_1.csv"),
        "pcs" => joinpath("data", "pcs.csv"),
        "call-threshold" => 0.8,  
        "%COMMAND%" => "from-actors", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => "final", 
        "batch-size" => nothing,
        "positivity-constraint" => 0.
    )
    bqtls = Symbol.(unique(CSV.read(parsed_args["from-actors"]["bqtls"], DataFrame).ID))
    tmle_inputs(parsed_args)

    ## Dataset file
    trait_data = CSV.read("final.data.csv", DataFrame)
    @test names(trait_data) == [
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2", "RSID_2", "RSID_102", 
        "RSID_17", "RSID_198", "RSID_99"]
    @test size(trait_data) == (490, 16)
    
    ## Parameter file: 
    outparameters = parameters_from_yaml("final.param.yaml")
    found_targets = Dict(
        :BINARY_1 => 0,
        :CONTINUOUS_2 => 0,
        :CONTINUOUS_1 => 0,
        :BINARY_2 => 0
    )
    for Ψ in outparameters
        if Ψ isa ATE
            ntreatments = length(Ψ.treatment)
            if ntreatments > 1
                @test all(Ψ.treatment[index].case == Ψ.treatment[index].control for index ∈ 2:ntreatments)
            end
        else
            @test Ψ isa IATE
            @test all(cc.case != cc.control for cc ∈ Ψ.treatment)
        end
        @test Ψ.covariates == [:COV_1, Symbol("21003"), Symbol("22001")]
        @test Ψ.confounders == [:PC1, :PC2]
        # The first treatment will be a bqtl
        @test keys(Ψ.treatment)[1] ∈ bqtls
        @test Ψ.treatment[1].case isa Int
        @test length(Ψ.treatment) ∈ [1, 2]
        found_targets[Ψ.target] += 1
    end
    # The number of parameters with various targets should be the same
    @test all(x == found_targets[:BINARY_1] for x in values(found_targets))
    # This is difficult to really check the ordering
    # Those correspond to the simple bQTL ATE
    first_treatments = keys(outparameters[1].treatment)
    @test all(keys(Ψ.treatment) == first_treatments for Ψ in outparameters[1:12])

    cleanup()
end

@testset "Test tmle_inputs from-actors: scenario 2" begin
    # Scenario:
    # - 2 sets of trans actors
    # - no extra treatment
    # - no extra covariate
    # - extra confounders
    # - orders 2 and 3
    # - no continuous phenotypes 
    # - batched
    parsed_args = Dict(
        "from-actors" => Dict{String, Any}(
            "bqtls" => joinpath("data", "bqtls.csv"), 
            "trans-actors-prefix" => joinpath("data", "trans_actors_2"),
            "extra-covariates" => nothing,
            "extra-treatments" => nothing,
            "extra-confounders" => joinpath("data", "extra_confounders.txt"),
            "orders" => "2,3",
            "genotypes-as-int" => false
            ),
        "traits" => joinpath("data", "traits_2.csv"),
        "pcs" => joinpath("data", "pcs.csv"),
        "call-threshold" => 0.8,  
        "%COMMAND%" => "from-actors", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => "final", 
        "batch-size" => 100,
        "positivity-constraint" => 0.,
    )
    bqtls = Symbol.(unique(CSV.read(parsed_args["from-actors"]["bqtls"], DataFrame).ID))
    
    tmle_inputs(parsed_args)
    
    ## Dataset file
    traits = CSV.read("final.data.csv", DataFrame)
    @test names(traits) == [
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "COV_1", "21003", "22001", 
        "PC1", "PC2", "RSID_2", "RSID_102", "RSID_17", "RSID_198", "RSID_99"]
    @test size(traits) == (490, 13)

    # Parameter files: 
    ouparameters_1 = parameters_from_yaml("final.param_1.yaml")
    @test size(ouparameters_1, 1) == 100
    ouparameters_2 = parameters_from_yaml("final.param_2.yaml")
    outparameters = vcat(ouparameters_1, ouparameters_2)
    
    found_targets = Dict(
        :BINARY_1 => 0,
        :BINARY_2 => 0
    )
    for Ψ in outparameters
        if Ψ isa ATE
            ntreatments = length(Ψ.treatment)
            @test all(Ψ.treatment[index].case == Ψ.treatment[index].control for index ∈ 2:ntreatments)
        else
            @test Ψ isa IATE
            @test all(cc.case != cc.control for cc ∈ Ψ.treatment)
        end
        @test Ψ.covariates == []
        @test Ψ.confounders == [:PC1, :PC2, :COV_1, Symbol("21003"), Symbol("22001")]
        # The first treatment will be a bqtl
        @test keys(Ψ.treatment)[1] ∈ bqtls
        @test Ψ.treatment[1].case isa String
        @test length(Ψ.treatment) ∈ [2, 3]
        found_targets[Ψ.target] += 1
    end
    # The number of parameters with various targets should be the same
    @test all(x == found_targets[:BINARY_1] for x in values(found_targets))
    # This is difficult to really check the ordering
    # Those correspond to the simple bQTL ATE
    first_treatments = keys(outparameters[1].treatment)
    @test all(keys(Ψ.treatment) == first_treatments for Ψ in outparameters[1:15])

    cleanup()

    # Adding positivity constraint, only 4 files are generated
    parsed_args["positivity-constraint"] = 0.1
    tmle_inputs(parsed_args)

    @test !isfile("final.param_2.yaml")
    outparameters = parameters_from_yaml("final.param_1.yaml")
    @test size(outparameters, 1) == 6
    cleanup()
end

end

true