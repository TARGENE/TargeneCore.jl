module TestFromActors

using Test
using CSV
using DataFrames
using TargeneCore
using YAML
using TMLE
using Arrow
using Serialization

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "tl_inputs", "test_utils.jl"))

function test_cates(Ψ, bqtls)
    nfixed = 0
    for T ∈ keys(Ψ.treatment_values)
        if T ∉ bqtls
            nfixed += 1
            @test Ψ.treatment_values[T].case == Ψ.treatment_values[T].control
        end
    end
    @test nfixed > 0
end

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
    settings = collect(TargeneCore.control_case_settings(TMLE.StatisticalATE, treatments, data))
    @test settings == [
        (T1 = (case = 0, control = 1),),
        (T1 = (case = 0, control = 2),),
        (T1 = (case = 1, control = 2),)
    ]
    # Two treatments IATEs
    treatments = [:T1, :T2]
    settings = collect(TargeneCore.control_case_settings(TMLE.StatisticalIATE, treatments, data))
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
    settings = collect(TargeneCore.control_case_settings(TMLE.StatisticalATE, treatments, data))
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
    settings = collect(TargeneCore.control_case_settings(TMLE.StatisticalIATE, treatments, data))
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
    settings = collect(TargeneCore.control_case_settings(TMLE.StatisticalATE, treatments, data))
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

    bqtl_file = joinpath(TESTDIR, "data", "bqtls_1.csv")
    trans_actors_prefix = joinpath(TESTDIR, "data", "trans_actors_1.csv")
    env_file = joinpath(TESTDIR, "data", "extra_treatments.txt")
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


@testset "Test tl_inputs from-actors: scenario 1" begin
    # Scenario:
    # - Trans-actors
    # - Extra Treatment
    # - Extra Covariates
    # - Order 1,2
    parsed_args = Dict(
        "out-prefix" => "final", 
        "batch-size" => nothing,
        "positivity-constraint" => 0.,
        "verbosity" => 0,

        "%COMMAND%" => "from-actors", 

        "from-actors" => Dict{String, Any}(
            "bqtls" => joinpath(TESTDIR, "data", "bqtls_1.csv"), 
            "trans-actors-prefix" => joinpath(TESTDIR, "data", "trans_actors_1.csv"),
            "extra-covariates" => joinpath(TESTDIR, "data", "extra_covariates.txt"),
            "extra-treatments" => joinpath(TESTDIR, "data", "extra_treatments.txt"),
            "extra-confounders" => nothing,
            "orders" => "1,2",
            "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
            "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
            "call-threshold" => 0.8,  
            "bgen-prefix" => joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb"), 
            ),
    )
    bqtls = Symbol.(unique(CSV.read(parsed_args["from-actors"]["bqtls"], DataFrame).ID))
    tl_inputs(parsed_args)

    ## Dataset file
    trait_data = DataFrame(Arrow.Table("final.data.arrow"))
    @test sort(names(trait_data)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2", "RSID_2", "RSID_102", 
        "RSID_17", "RSID_198", "RSID_99"])
    @test size(trait_data) == (490, 16)
    
    ## Output estimands: 
    output_estimands = deserialize("final.estimands_1.jls").estimands
    found_targets = Dict(
        :BINARY_1 => 0,
        :CONTINUOUS_2 => 0,
        :CONTINUOUS_1 => 0,
        :BINARY_2 => 0
    )
    for Ψ in output_estimands
        if Ψ isa TMLE.StatisticalATE
            # These are CATEs
            if length(Ψ.treatment_values) > 1
                test_cates(Ψ, bqtls)
            end
        else
            @test Ψ isa TMLE.StatisticalIATE
            @test all(cc.case != cc.control for cc ∈ Ψ.treatment_values)
        end
        @test Ψ.outcome_extra_covariates == (Symbol("21003"), Symbol("22001"), :COV_1)
        @test all(x == (:PC1, :PC2) for x ∈ Ψ.treatment_confounders)
        @test any(T ∈ bqtls for T ∈ keys(Ψ.treatment_values))
        @test length(Ψ.treatment_values) ∈ [1, 2]
        found_targets[Ψ.outcome] += 1
    end
    # The number of estimands with various targets should be the same
    @test all(x == found_targets[:BINARY_1] for x in values(found_targets))
    # This is difficult to really check the full ordering
    # Only check the first estimands share the same propensity score
    first_treatments = keys(output_estimands[1].treatment_values)
    @test all(keys(Ψ.treatment_values) == first_treatments for Ψ in output_estimands[1:12])

    cleanup()
end

@testset "Test tl_inputs from-actors: scenario 2" begin
    # Scenario:
    # - 2 sets of trans actors
    # - no extra treatment
    # - no extra covariate
    # - extra confounders
    # - orders 2 and 3
    # - no continuous phenotypes 
    # - batched
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => "final", 
        "batch-size" => 100,
        "positivity-constraint" => 0.,

        "%COMMAND%" => "from-actors", 

        "from-actors" => Dict{String, Any}(
            "bqtls" => joinpath(TESTDIR, "data", "bqtls_1.csv"), 
            "trans-actors-prefix" => joinpath(TESTDIR, "data", "trans_actors_2"),
            "extra-covariates" => nothing,
            "extra-treatments" => nothing,
            "extra-confounders" => joinpath(TESTDIR, "data", "extra_confounders.txt"),
            "orders" => "2,3",
            "traits" => joinpath(TESTDIR, "data", "traits_2.csv"),
            "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
            "call-threshold" => 0.8,  
            "bgen-prefix" => joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb"), 
        ),
    )
    bqtls = Symbol.(unique(CSV.read(parsed_args["from-actors"]["bqtls"], DataFrame).ID))
    
    tl_inputs(parsed_args)
    
    ## Dataset file
    traits = DataFrame(Arrow.Table("final.data.arrow"))
    @test sort(names(traits)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "COV_1", "21003", "22001", 
        "PC1", "PC2", "RSID_2", "RSID_102", "RSID_17", "RSID_198", "RSID_99"])
    @test size(traits) == (490, 13)

    # Parameter files: 
    output_estimands_1 = deserialize("final.estimands_1.jls").estimands
    @test size(output_estimands_1, 1) == 100
    output_estimands_2 = deserialize("final.estimands_2.jls").estimands
    output_estimands = vcat(output_estimands_1, output_estimands_2)
    
    found_targets = Dict(
        :BINARY_1 => 0,
        :BINARY_2 => 0
    )
    for Ψ in output_estimands
        if Ψ isa TMLE.StatisticalATE
            nfixed = 0
            for T ∈ keys(Ψ.treatment_values)
                if T ∉ bqtls
                    nfixed += 1
                    @test Ψ.treatment_values[T].case == Ψ.treatment_values[T].control
                end
            end
            @test nfixed > 0
        else
            @test Ψ isa TMLE.StatisticalIATE
            @test all(cc.case != cc.control for cc ∈ Ψ.treatment_values)
        end
        @test Ψ.outcome_extra_covariates == ()
        @test all(x == (Symbol("21003"), Symbol("22001"), :COV_1, :PC1, :PC2) for x ∈ Ψ.treatment_confounders)
        @test any(T ∈ bqtls for T ∈ keys(Ψ.treatment_values))
        @test length(Ψ.treatment_values) ∈ [2, 3]
        found_targets[Ψ.outcome] += 1
    end
    # The number of estimands with various targets should be the same
    @test all(x == found_targets[:BINARY_1] for x in values(found_targets))
    # This is difficult to really check the ordering
    first_treatments = keys(output_estimands[1].treatment_values)
    @test all(keys(Ψ.treatment_values) == first_treatments for Ψ in output_estimands[1:15])

    cleanup()

    # Adding positivity constraint, only 4 files are generated
    parsed_args["positivity-constraint"] = 0.1
    tl_inputs(parsed_args)

    @test !isfile("final.param_2.yaml")
    output_estimands = deserialize("final.estimands_1.jls").estimands
    @test size(output_estimands, 1) == 6
    cleanup()
end

@testset "Test tl_inputs from-actors: scenario 3" begin
    # Scenario:
    # - Trans-actors
    # - Extra Treatment
    # - Extra Covariates
    # - Order 1,2
    # - More than 1 TF present
    parsed_args = Dict(
        "verbosity" => 0,
        "out-prefix" => "final", 
        "batch-size" => nothing,
        "positivity-constraint" => 0.,

        "%COMMAND%" => "from-actors", 

        "from-actors" => Dict{String, Any}(
            "bqtls" => joinpath(TESTDIR, "data", "bqtls_2.csv"), 
            "trans-actors-prefix" => joinpath(TESTDIR, "data", "trans_actors_3.csv"),
            "extra-covariates" => joinpath(TESTDIR, "data", "extra_covariates.txt"),
            "extra-treatments" => joinpath(TESTDIR, "data", "extra_treatments.txt"),
            "extra-confounders" => nothing,
            "orders" => "1,2",
            "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
            "pcs" => joinpath(TESTDIR, "data", "pcs.csv"),
            "call-threshold" => 0.8,  
            "bgen-prefix" => joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb"),
            ), 
    )
    bqtls = Symbol.(unique(CSV.read(parsed_args["from-actors"]["bqtls"], DataFrame).ID))
    tl_inputs(parsed_args)

    ## Dataset file
    trait_data = DataFrame(Arrow.Table("final.data.arrow"))
    @test sort(names(trait_data)) == sort([
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", 
        "COV_1", "21003", "22001", "TREAT_1", "PC1", "PC2", "RSID_2", "RSID_102", 
        "RSID_17", "RSID_198", "RSID_99"])
    @test size(trait_data) == (490, 16)
    
    ## Parameter file: 
    tf_output_estimands = [
        deserialize("final.TF1.estimands_1.jls").estimands, 
        deserialize("final.TF2.estimands_1.jls").estimands
    ]
    found_targets = Dict(
        :BINARY_1 => 0,
        :CONTINUOUS_2 => 0,
        :CONTINUOUS_1 => 0,
        :BINARY_2 => 0
    )
    tf_transactors = [[:RSID_102], [:RSID_2]]
    tf_bqtls = [[:RSID_17, :RSID_99], [:RSID_17, :RSID_198]]
    bqtl_transactor_count = 0
    for tf in (1,2)
        for Ψ in tf_output_estimands[tf]
            @test any(bqtl ∈ keys(Ψ.treatment_values) for bqtl ∈ tf_bqtls[tf])
            if length(Ψ.treatment_values) > 1
                if any(transactor ∈ keys(Ψ.treatment_values) for transactor ∈ tf_transactors[tf])
                    bqtl_transactor_count += 1
                end
            end
        end
    end
    @test bqtl_transactor_count > 100

    cleanup()
end

end

true