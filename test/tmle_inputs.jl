module TestTMLEInputs

using Test
using CSV
using DataFrames
using TargeneCore
using YAML

function cleanup()
    for file in readdir()
        if startswith(file, "final.")
            rm(file)
        end
    end
end

@testset "Test tmle_inputs with-param-files: scenario 1" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - genetic and extra confounders
    # - no covariates
    # - no extra treatments
    # - no batch size
    parsed_args = Dict(
        "with-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_")), 
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => nothing, 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => joinpath("data", "extra_confounders.csv"), 
        "%COMMAND%" => "with-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => nothing,
        "phenotype-batch-size" => nothing,
        "positivity-constraint" => 0.01
    )
    @test_logs(
        (:warn, "Some treatment variables could not be read from the data files and associated parameter files will not be processed: TREAT_1"), 
        tmle_inputs(parsed_args)
    )

    # Data Files
    confounders = CSV.read("final.confounders.csv", DataFrame)
    @test names(confounders) == ["SAMPLE_ID", "PC1", "PC2", "21003", "22001"]
    @test size(confounders) == (490, 5)

    treatments = CSV.read("final.treatments.csv", DataFrame)
    @test size(treatments) == (490, 3)
    @test names(treatments) == ["SAMPLE_ID", "RSID_2", "RSID_198"]
    @test Set(unique(treatments[!, "RSID_2"])) == Set([missing, 0, 1, 2])
    @test Set(unique(treatments[!, "RSID_198"])) == Set([0, 1, 2])

    binary_phenotypes = CSV.read("final.binary-phenotypes.csv", DataFrame)
    @test size(binary_phenotypes) == (490, 3)
    @test names(binary_phenotypes) == ["SAMPLE_ID", "BINARY_1", "BINARY_2"]

    continuous_phenotypes = CSV.read("final.continuous-phenotypes.csv", DataFrame)  
    @test size(continuous_phenotypes) == (490, 3)
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "CONTINUOUS_1", "CONTINUOUS_2"]

    # Parameter files: untouched because no phenotype batches were specified
    # They are deduplicated for both continuous and binary phenotypes
    @test YAML.load_file(joinpath("config", "param_1.yaml")) == YAML.load_file("final.binary.parameter_1.yaml") == YAML.load_file("final.continuous.parameter_1.yaml")
    @test YAML.load_file(joinpath("config", "param_2.yaml")) == YAML.load_file("final.binary.parameter_2.yaml") == YAML.load_file("final.continuous.parameter_2.yaml")
    
    cleanup()
end

@testset "Test tmle_inputs with-param-files: scenario 2" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - no extra confounders
    # - covariates
    # - SNP and extra treatments
    # - batch size: 1
    parsed_args = Dict(
        "with-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_1")), 
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => joinpath("data", "extra_treatments.csv"), 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => nothing, 
        "%COMMAND%" => "with-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => joinpath("data", "covariates.csv"),
        "phenotype-batch-size" => 1,
        "positivity-constraint" => 0.01
    )
    tmle_inputs(parsed_args)

    confounders = CSV.read("final.confounders.csv", DataFrame)
    @test names(confounders) == ["SAMPLE_ID", "PC1", "PC2"]
    @test size(confounders) == (490, 3)

    treatments = CSV.read("final.treatments.csv", DataFrame)
    @test size(treatments) == (490, 3)
    @test names(treatments) == ["SAMPLE_ID", "RSID_2", "TREAT_1"]

    binary_phenotypes = CSV.read("final.binary-phenotypes.csv", DataFrame)
    @test size(binary_phenotypes) == (490, 3)
    @test names(binary_phenotypes) == ["SAMPLE_ID", "BINARY_1", "BINARY_2"]

    continuous_phenotypes = CSV.read("final.continuous-phenotypes.csv", DataFrame)  
    @test size(continuous_phenotypes) == (490, 3)
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "CONTINUOUS_1", "CONTINUOUS_2"]

    covariates = CSV.read("final.covariates.csv", DataFrame)
    @test names(covariates) == ["SAMPLE_ID", "COV_1"]
    @test size(covariates) == (490, 2)
    # Parameter files: modified because phenotype batches was specified
    # They are deduplicated for both continuous and binary phenotypes
    origin_1 = YAML.load_file(joinpath("config", "param_1.yaml"))
    binary_1_1 = YAML.load_file("final.binary.parameter_1.yaml")
    binary_1_2 = YAML.load_file("final.binary.parameter_2.yaml")
    continuous_1_1 = YAML.load_file("final.continuous.parameter_1.yaml")
    continuous_1_2 = YAML.load_file("final.continuous.parameter_2.yaml")
    # Parameters and Treatments sections unchanged
    @test binary_1_1["Parameters"] == binary_1_2["Parameters"] == continuous_1_1["Parameters"] == origin_1["Parameters"]
    @test binary_1_1["Treatments"] == binary_1_2["Treatments"] == continuous_1_2["Treatments"] == origin_1["Treatments"]
    # Phenotypes sections changed
    @test binary_1_1["Phenotypes"] == ["BINARY_1"]
    @test binary_1_2["Phenotypes"] == ["BINARY_2"]
    @test continuous_1_1["Phenotypes"] == ["CONTINUOUS_1"]
    @test continuous_1_2["Phenotypes"] == ["CONTINUOUS_2"]

    origin_2 = YAML.load_file(joinpath("config", "param_1_with_extra_treatment.yaml"))
    binary_2_1 = YAML.load_file("final.binary.parameter_3.yaml")
    binary_2_2 = YAML.load_file("final.binary.parameter_4.yaml")
    continuous_2_1 = YAML.load_file("final.continuous.parameter_3.yaml")
    continuous_2_2 = YAML.load_file("final.continuous.parameter_4.yaml")
    # Parameters and Treatments sections unchanged
    @test binary_2_1["Parameters"] == binary_2_2["Parameters"] == continuous_2_1["Parameters"] == origin_2["Parameters"]
    @test binary_2_1["Treatments"] == binary_2_2["Treatments"] == continuous_2_2["Treatments"] == origin_2["Treatments"]
    # Phenotypes sections changed
    @test binary_2_1["Phenotypes"] == ["BINARY_1"]
    @test binary_2_2["Phenotypes"] == ["BINARY_2"]
    @test continuous_2_1["Phenotypes"] == ["CONTINUOUS_1"]
    @test continuous_2_2["Phenotypes"] == ["CONTINUOUS_2"]

    cleanup()
end

@testset "Test tmle_inputs with-asb-trans: scenario 1" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - genetic and extra confounders
    # - no covariates
    # - no extra treatments
    # - no batch size
    # - no param
    parsed_args = Dict(
        "with-asb-trans" => Dict{String, Any}(
            "asb-prefix" => joinpath("data", "asb_files", "asb"), 
            "trans-actors" => joinpath("data", "trans_actors_fake.csv"),
            "param-prefix" => nothing
            ),
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => nothing, 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => joinpath("data", "extra_confounders.csv"), 
        "%COMMAND%" => "with-asb-trans", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => nothing,
        "phenotype-batch-size" => nothing,
        "positivity-constraint" => 0.01
    )
    tmle_inputs(parsed_args)

    cleanup()
end

end

true