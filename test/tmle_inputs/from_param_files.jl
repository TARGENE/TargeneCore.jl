module TestFromParamFiles

using Test
using CSV
using DataFrames
using TargeneCore
using YAML


include("test_utils.jl")


#####################################################################
###############               UNIT TESTS              ###############
#####################################################################

@testset "Test snps_and_variables_from_param_files" begin
    pcs = TargeneCore.read_data(joinpath("data", "pcs.csv"))
    traits = TargeneCore.read_data(joinpath("data", "traits_1.csv"))
    # extraW, extraT, extraC are parsed from all param_files
    param_files = TargeneCore.load_param_files(joinpath("config", "param"))
    snps, variables = TargeneCore.snps_and_variables_from_param_files(param_files, pcs, traits)
    @test snps == ["RSID_2", "RSID_198"]
    @test variables == Dict(
        "Y"           => ["BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2"],
        "extraW"      => ["22001"],
        "extraT"      => ["TREAT_1"],
        "PCs"         => ["PC1", "PC2"],
        "extraC"      => ["COV_1", "21003", "22001"]
    )

    # In the following scenario, some variables are interpreted as targets while they shouldnt
    # because no parameter file describes them otherwise. This does not lead to any error since
    # those variables can be interpreted as either binary or continuous variables.
    param_files = TargeneCore.load_param_files(joinpath("config", "param_1_with"))
    snps, variables = TargeneCore.snps_and_variables_from_param_files(param_files, pcs, traits)
    @test snps == ["RSID_2"]
    @test variables == Dict(
        "Y"           => ["BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", "COV_1", "21003", "22001",],
        "extraW"      => [],
        "extraT"      => ["TREAT_1"],
        "PCs"         => ["PC1", "PC2"],
        "extraC"      => [],
    )

end

#####################################################################
###############           END-TO-END TESTS            ###############
#####################################################################


@testset "Test tmle_inputs from-param-files: scenario 1" begin
    parsed_args = Dict(
        "from-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_")), 
        "traits" => joinpath("data", "traits_1.csv"),
        "pcs" => joinpath("data", "pcs.csv"),
        "call-threshold" => 0.8, 
        "%COMMAND%" => "from-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => "final", 
        "phenotype-batch-size" => nothing,
        "positivity-constraint" => 0.
    )

    tmle_inputs(parsed_args)

    # Data Files
    data = CSV.read("final.data.csv", DataFrame)
    @test names(data) == [
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", 
        "CONTINUOUS_2", "COV_1", "21003", "22001", "TREAT_1", 
        "PC1", "PC2", "RSID_2", "RSID_198"
    ]
    # Parameter files:
    # - 2 of them specify the targets
    # - The rest should incorporate all targets
    for index in 1:5
        param_file = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], index))
        if param_file["T"] == ["RSID_2"] && param_file["W"] == ["PC1", "PC2"]
            if haskey(param_file, "C")
                @test param_file["C"] == ["COV_1", 21003]
                @test param_file["Y"] == ["BINARY_1"]
            else
                @test param_file["Y"] == ["BINARY_1", "CONTINUOUS_1"]
            end
        else
            @test param_file["Y"] == ["BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2"]
            if param_file["T"] == ["RSID_2", "TREAT_1"]
                @test param_file["W"] == ["PC1", "PC2"]
                @test !haskey(param_file, "C")
            elseif param_file["T"] == ["RSID_2"]
                @test param_file["W"] == [22001, "PC1", "PC2"]
                @test param_file["C"] == ["COV_1", 21003]
            else
                @test param_file["T"] == ["RSID_2", "RSID_198"]
                @test param_file["W"] == ["PC1", "PC2"]
                @test param_file["C"] == [22001]
            end
        end
    end

    cleanup()
end

@testset "Test tmle_inputs from-param-files: scenario 2" begin
    parsed_args = Dict(
        "from-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_2")), 
        "traits" => joinpath("data", "traits_1.csv"),
        "pcs" => joinpath("data", "pcs.csv"),
        "call-threshold" => 0.8, 
        "%COMMAND%" => "from-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => "final", 
        "phenotype-batch-size" => nothing,
        "positivity-constraint" => 0.
    )

    tmle_inputs(parsed_args)
    param_file = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 1))
    @test param_file == Dict(
        "Y"          => ["BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", "COV_1", "21003", "TREAT_1"],
        "W"          => ["PC1", "PC2"],
        "T"          => ["RSID_2", "RSID_198"],
        "C"          => [22001],
        "Parameters" => [Dict("name"=>"ATE", "RSID_2"=>Dict("case"=>1, "control"=>0), "RSID_198"=>Dict("case"=>1, "control"=>0))]
    )

    cleanup()
end

@testset "Test tmle_inputs from-param-files: scenario 3" begin
    # Batched
    parsed_args = Dict(
        "from-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_1_with")), 
        "traits" => joinpath("data", "traits_2.csv"),
        "pcs" => joinpath("data", "pcs.csv"),
        "call-threshold" => 0.8, 
        "%COMMAND%" => "from-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => "final", 
        "phenotype-batch-size" => 2,
        "positivity-constraint" => 0.
    )
    tmle_inputs(parsed_args)
    # 2 first variables in a batch
    param_file_1 = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 1))
    @test param_file_1 == Dict(
        "Y"          => ["BINARY_1", "BINARY_2"],
        "W"          => ["PC1", "PC2"],
        "T"          => ["RSID_2", "TREAT_1"],
        "Parameters" => Dict[Dict("name"=>"IATE", "TREAT_1"=>Dict("case"=>1, "control"=>0), "RSID_2"=>Dict("case"=>1, "control"=>0))]
    )
    # 22001 is interpreted as a binary target
    param_file_2 = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 2))
    @test param_file_2 == Dict(
        "Y"          => ["COV_1", "21003"],
        "W"          => ["PC1", "PC2"],
        "T"          => ["RSID_2", "TREAT_1"],
        "Parameters" => Dict[Dict("name"=>"IATE", "TREAT_1"=>Dict("case"=>1, "control"=>0), "RSID_2"=>Dict("case"=>1, "control"=>0))]
    )
    # COV_1 & 21003 are interpreted as continuous targets
    param_file_3 = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 3))
    @test param_file_3 == Dict(
        "Y"          => ["22001"],
        "W"          => ["PC1", "PC2"],
        "T"          => ["RSID_2", "TREAT_1"],
        "Parameters" => Dict[Dict("name"=>"IATE", "TREAT_1"=>Dict("case"=>1, "control"=>0), "RSID_2"=>Dict("case"=>1, "control"=>0))]
    )

    cleanup()
end

@testset "Test tmle_inputs from-param-files: scenario 4" begin
    parsed_args = Dict(
        "from-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_with_specified_targets")), 
        "traits" => joinpath("data", "traits_1.csv"),
        "pcs" => joinpath("data", "pcs.csv"),
        "call-threshold" => 0.8, 
        "%COMMAND%" => "from-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "out-prefix" => "final", 
        "phenotype-batch-size" => nothing,
        "positivity-constraint" => 0.
    )
    tmle_inputs(parsed_args)

    # One parameter file has no extra covariates
    param_file_1 = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 1))
    @test param_file_1 == Dict(
        "Y"          => ["BINARY_1", "CONTINUOUS_1"],
        "W"          => ["PC1", "PC2"],
        "T"          => ["RSID_2"],
        "Parameters" => [Dict("name"=>"ATE", "RSID_2"=>Dict("case"=>1, "control"=>0))]
    )

    # originating from a parameter file with extra covariates
    param_file_2 = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 2))
    @test param_file_2 == Dict(
        "Y"          => ["BINARY_1"],
        "W"          => ["PC1", "PC2"],
        "T"          => ["RSID_2"],
        "C"          => ["COV_1", 21003],
        "Parameters" => [Dict("name"=>"ATE", "RSID_2"=>Dict("case"=>1, "control"=>0))]
    )
    cleanup()
end

end

true