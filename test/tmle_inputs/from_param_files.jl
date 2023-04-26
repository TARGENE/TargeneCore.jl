module TestFromParamFiles

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

@testset "Test get_variables" begin
    traits = TargeneCore.read_data(joinpath("data", "traits_1.csv"))
    pcs = TargeneCore.read_data(joinpath("data", "pcs.csv"))
    # extraW, extraT, extraC are parsed from all param_files
    param_groups = parameters_from_yaml(joinpath("config", "parameters_integer_encoding.yaml"))
    variables = TargeneCore.get_variables(param_groups, traits, pcs)
    @test variables.variants == Set([:RSID_198, :RSID_2])
    @test variables.targets == Set([:BINARY_1, :CONTINUOUS_2, :CONTINUOUS_1, :BINARY_2])
    @test variables.pcs == Set([:PC1, :PC2])
end

@testset "Test get_genotype_encoding" begin
    variants = Set([:RSID_2, :RSID_198])
    # Throwing because one is Integer representation and the other a String
    parameters = parameters_from_yaml(joinpath("config", "parameters_corrupted_encodings.yaml"))
    @test_throws TargeneCore.MismatchedCaseControlEncodingError() TargeneCore.get_genotype_encoding(parameters, variants)
    # Only string is fine
    parameters = parameters_from_yaml(joinpath("config", "parameters_string_encoding.yaml"))
    genotypes_asint = TargeneCore.get_genotype_encoding(parameters, variants)
    @test genotypes_asint == false
    # Only Int is fine
    parameters = parameters_from_yaml(joinpath("config", "parameters_integer_encoding.yaml"))
    genotypes_asint = TargeneCore.get_genotype_encoding(parameters, variants)
    @test genotypes_asint == true
end

@testset "Test adjust_parameter_sections" begin
    ## With string encoded variants
    genotypes = DataFrame(
        SAMPLE_ID = [1, 2, 3],
        RSID_198 = ["GA", "GG", "AA"],
        RSID_2 = ["AG", "GG", "AA"],
    )
    pcs = Set([:PC1, :PC2])
    variants_alleles = Dict(:RSID_198 => Set(genotypes.RSID_198))
    # AG is not in the genotypes bu GA is
    Ψ = parameters_from_yaml(joinpath("config", "parameters_string_encoding.yaml"))[4]
    @test Ψ.treatment.RSID_198 == (case="AG", control="AA")
    new_Ψ = TargeneCore.adjust_parameter_sections(Ψ, variants_alleles, pcs)
    @test new_Ψ.target == Ψ.target
    @test new_Ψ.covariates == Ψ.covariates
    @test new_Ψ.confounders == [:PC1, :PC2]
    @test new_Ψ.treatment == (
        RSID_2 = (case = "AA", control = "GG"),
        RSID_198 = (case = "GA", control = "AA")
    )

    # If the allele is not present 
    variants_alleles = Dict(:RSID_198 => Set(["AA"]))
    @test_throws TargeneCore.AbsentAlleleError("RSID_198", "AG") TargeneCore.adjust_parameter_sections(Ψ, variants_alleles, pcs)
    
    ## With integer encoded variants
    genotypes = DataFrame(
        SAMPLE_ID = [1, 2, 3],
        RSID_2 = [1, 2, 0],
    )
    variants_alleles = Dict(:RSID_2 => Set(genotypes.RSID_2))
    Ψ = TargeneCore.parameters_from_yaml(joinpath("config", "parameters_integer_encoding.yaml"))[1]
    new_Ψ = TargeneCore.adjust_parameter_sections(Ψ, variants_alleles, pcs)
    @test new_Ψ.target == Ψ.target
    @test new_Ψ.covariates == Ψ.covariates
    @test new_Ψ.confounders == [:PC1, :PC2]
    @test new_Ψ.treatment == Ψ.treatment

    variants_alleles = Dict(:RSID_2 => Set([]))
    @test_throws TargeneCore.AbsentAlleleError("RSID_2", 1) TargeneCore.adjust_parameter_sections(Ψ, variants_alleles, pcs)
end

#####################################################################
###############           END-TO-END TESTS            ###############
#####################################################################


@testset "Test tmle_inputs from-param-files: parameters_string_encoding.yaml" begin
    # Genotypes encoded as strings
    # No batching of parameter files
    # No positivity constraint
    parsed_args = Dict(
        "from-param-files" => Dict{String, Any}("paramfile" => joinpath("config", "parameters_string_encoding.yaml")), 
        "traits" => joinpath("data", "traits_1.csv"),
        "pcs" => joinpath("data", "pcs.csv"),
        "call-threshold" => 0.8, 
        "%COMMAND%" => "from-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => "final", 
        "phenotype-batch-size" => nothing,
        "positivity-constraint" => 0.,
    )

    tmle_inputs(parsed_args)

    ## Data File
    data = CSV.read("final.data.csv", DataFrame)
    @test names(data) == [
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", 
        "CONTINUOUS_2", "COV_1", "21003", "22001", "TREAT_1", 
        "PC1", "PC2", "RSID_2", "RSID_198"
    ]
    @test size(data) == (490, 13)

    ## Parameter file:
    # There are 5 initial parameters containing a *
    # Those are duplicated for each of the 4 targets.
    # The PCs are appended to the confounders.
    # GA is corrected to AG to match the data
    outparameters = parameters_from_yaml("final.param.yaml")
    # Input parameters 2 and 3 (ATE and CM) share treatment, confounders, covariates and will be grouped together
    # Input parameters 4 and 5 (ATE and CM again) share treatment, confounders, covariates and will be grouped together
    expected_parameters = [
        ATE(:BINARY_1, (RSID_2 = (case = "AA", control = "GG"),), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
        CM(:BINARY_1, (RSID_2 = "AA",), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
        ATE(:BINARY_2, (RSID_2 = (case = "AA", control = "GG"),), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
        CM(:BINARY_2, (RSID_2 = "AA",), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
        ATE(:CONTINUOUS_1, (RSID_2 = (case = "AA", control = "GG"),), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
        CM(:CONTINUOUS_1, (RSID_2 = "AA",), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
        ATE(:CONTINUOUS_2, (RSID_2 = (case = "AA", control = "GG"),), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
        CM(:CONTINUOUS_2, (RSID_2 = "AA",), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
        ATE(:BINARY_1, (RSID_2 = (case = "AA", control = "GG"), RSID_198 = (case = "AG", control = "AA")), [:PC1, :PC2], [Symbol("22001")]),
        CM(:BINARY_1, (RSID_2 = "GG", RSID_198 = "AG"), [:PC1, :PC2], [Symbol("22001")]),
        ATE(:BINARY_2, (RSID_2 = (case = "AA", control = "GG"), RSID_198 = (case = "AG", control = "AA")), [:PC1, :PC2], [Symbol("22001")]),
        CM(:BINARY_2, (RSID_2 = "GG", RSID_198 = "AG"), [:PC1, :PC2], [Symbol("22001")]),
        ATE(:CONTINUOUS_1, (RSID_2 = (case = "AA", control = "GG"), RSID_198 = (case = "AG", control = "AA")), [:PC1, :PC2], [Symbol("22001")]),
        CM(:CONTINUOUS_1, (RSID_2 = "GG", RSID_198 = "AG"), [:PC1, :PC2], [Symbol("22001")]),
        ATE(:CONTINUOUS_2, (RSID_2 = (case = "AA", control = "GG"), RSID_198 = (case = "AG", control = "AA")), [:PC1, :PC2], [Symbol("22001")]),
        CM(:CONTINUOUS_2, (RSID_2 = "GG", RSID_198 = "AG"), [:PC1, :PC2], [Symbol("22001")]),
        IATE(:BINARY_1, (RSID_2 = (case = "AA", control = "GG"), TREAT_1 = (case = 1, control = 0)), [:PC1, :PC2], Symbol[]),
        IATE(:BINARY_2, (RSID_2 = (case = "AA", control = "GG"), TREAT_1 = (case = 1, control = 0)), [:PC1, :PC2], Symbol[]),
        IATE(:CONTINUOUS_1, (RSID_2 = (case = "AA", control = "GG"), TREAT_1 = (case = 1, control = 0)), [:PC1, :PC2], Symbol[]),
        IATE(:CONTINUOUS_2, (RSID_2 = (case = "AA", control = "GG"), TREAT_1 = (case = 1, control = 0)), [:PC1, :PC2], Symbol[]),
    ]
    @test outparameters == expected_parameters
    cleanup()

    # Increase positivity constraint
    parsed_args["positivity-constraint"] = 0.01
    tmle_inputs(parsed_args)
    # The IATES are the most sensitives
    outparameters = parameters_from_yaml("final.param.yaml")
    @test all(Ψ isa Union{CM, ATE} for Ψ in outparameters)
    @test size(outparameters, 1) == 16

    cleanup()

    parsed_args["positivity-constraint"] = 1.
    @test_throws TargeneCore.NoRemainingParamsError(1.) tmle_inputs(parsed_args)

end

@testset "Test tmle_inputs from-param-files: parameters_integer_encoding.yaml" begin
    parsed_args = Dict(
        "from-param-files" => Dict{String, Any}("paramfile" => joinpath("config", "parameters_integer_encoding.yaml")), 
        "traits" => joinpath("data", "traits_1.csv"),
        "pcs" => joinpath("data", "pcs.csv"),
        "call-threshold" => 0.8, 
        "%COMMAND%" => "from-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => "final", 
        "phenotype-batch-size" => 8,
        "positivity-constraint" => 0.,
    )

    tmle_inputs(parsed_args)
    ## Data File
    data = CSV.read("final.data.csv", DataFrame)
    @test names(data) == [
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", 
        "CONTINUOUS_2", "COV_1", "21003", "22001", "TREAT_1", 
        "PC1", "PC2", "RSID_2", "RSID_198"
    ]
    @test size(data) == (490, 13)

    ## Parameter files:
    # There are 5 initial parameters containing a *
    # Those are duplicated for each of the 4 targets.
    # The PCs are appended to the confounders.
    # Parameters are batched by 6
    expected_parameters = Dict(
        1 => [
            ATE(:BINARY_1, (RSID_2 = (case = 1, control = 0),), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
            CM(:BINARY_1, (RSID_2 = 1,), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
            ATE(:BINARY_2, (RSID_2 = (case = 1, control = 0),), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
            CM(:BINARY_2, (RSID_2 = 1,), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
            ATE(:CONTINUOUS_1, (RSID_2 = (case = 1, control = 0),), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
            CM(:CONTINUOUS_1, (RSID_2 = 1,), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
            ATE(:CONTINUOUS_2, (RSID_2 = (case = 1, control = 0),), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
            CM(:CONTINUOUS_2, (RSID_2 = 1,), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")])
        ],
        2 => [
            ATE(:BINARY_1, (RSID_2 = (case = 1, control = 0), RSID_198 = (case = 1, control = 0)), [:PC1, :PC2], [Symbol("22001")]),
            CM(:BINARY_1, (RSID_2 = 0, RSID_198 = 0), [:PC1, :PC2], [Symbol("22001")]),
            ATE(:BINARY_2, (RSID_2 = (case = 1, control = 0), RSID_198 = (case = 1, control = 0)), [:PC1, :PC2], [Symbol("22001")]),
            CM(:BINARY_2, (RSID_2 = 0, RSID_198 = 0), [:PC1, :PC2], [Symbol("22001")]),
            ATE(:CONTINUOUS_1, (RSID_2 = (case = 1, control = 0), RSID_198 = (case = 1, control = 0)), [:PC1, :PC2], [Symbol("22001")]),
            CM(:CONTINUOUS_1, (RSID_2 = 0, RSID_198 = 0), [:PC1, :PC2], [Symbol("22001")]),
            ATE(:CONTINUOUS_2, (RSID_2 = (case = 1, control = 0), RSID_198 = (case = 1, control = 0)), [:PC1, :PC2], [Symbol("22001")]),
            CM(:CONTINUOUS_2, (RSID_2 = 0, RSID_198 = 0), [:PC1, :PC2], [Symbol("22001")])
        ],
        3 => [
            IATE(:BINARY_1, (RSID_2 = (case = 1, control = 0), TREAT_1 = (case = 1, control = 0)), [:PC1, :PC2], Symbol[]),
            IATE(:BINARY_2, (RSID_2 = (case = 1, control = 0), TREAT_1 = (case = 1, control = 0)), [:PC1, :PC2], Symbol[]),
            IATE(:CONTINUOUS_1, (RSID_2 = (case = 1, control = 0), TREAT_1 = (case = 1, control = 0)), [:PC1, :PC2], Symbol[]),
            IATE(:CONTINUOUS_2, (RSID_2 = (case = 1, control = 0), TREAT_1 = (case = 1, control = 0)), [:PC1, :PC2], Symbol[])
        ]
    )
    for index in 1:3
        outparameters = parameters_from_yaml(string("final.param_", index, ".yaml"))
        @test outparameters == expected_parameters[index]
    end
    cleanup()
end

@testset "Test tmle_inputs from-param-files: param_without_wildcards.yaml" begin
    parsed_args = Dict(
        "from-param-files" => Dict{String, Any}("paramfile" => joinpath("config", "param_without_wildcards.yaml")), 
        "traits" => joinpath("data", "traits_1.csv"),
        "pcs" => joinpath("data", "pcs.csv"),
        "call-threshold" => 0.8, 
        "%COMMAND%" => "from-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => "final", 
        "phenotype-batch-size" => 2,
        "positivity-constraint" => 0.,
    )
    tmle_inputs(parsed_args)
    
    ## Data File
    data = CSV.read("final.data.csv", DataFrame)
    @test names(data) == [
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", 
        "CONTINUOUS_2", "COV_1", "21003", "22001", "TREAT_1", 
        "PC1", "PC2", "RSID_2"
    ]
    @test size(data) == (490, 12)

    ## Parameter files:
    # There are 3 initial parameters, 1 containing a *
    # that will be duplicated for each of the 4 targets.
    # The PCs are appended to the confounders.
    # Parameters are batched by 2
    expected_parameters = Dict(
        1 => [
            CM(:BINARY_1, (RSID_2 = 1,), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
            CM(:BINARY_2, (RSID_2 = 1,), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")])
        ],
        2 => [
            CM(:CONTINUOUS_1, (RSID_2 = 1,), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
            ATE(:CONTINUOUS_2, (RSID_2 = (case = 1, control = 0),), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")])
        ],
        3 => [
            CM(:CONTINUOUS_2, (RSID_2 = 1,), [Symbol("22001"), :PC1, :PC2], [:COV_1, Symbol("21003")]),
            IATE(:BINARY_1, (RSID_2 = (case = 1, control = 0), TREAT_1 = (case = 1, control = 0)), [:PC1, :PC2], Symbol[])
        ]
    )
    for index in 1:3
        outparameters = parameters_from_yaml(string("final.param_", index, ".yaml"))
        @test outparameters == expected_parameters[index]
    end

    cleanup()
end


end

true