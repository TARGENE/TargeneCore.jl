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

@testset "Test get_genotype_encoding" begin
    # Throwing because one is Integer representation and the other a String
    param_files = TargeneCore.load_param_files(joinpath("config", "param_2"))
    @test_throws TargeneCore.MismatchedCaseControlEncodingError() TargeneCore.get_genotype_encoding(param_files, ["RSID_2", "RSID_198"])
    # Only string is fine
    param_files = TargeneCore.load_param_files(joinpath("config", "param_2_with"))
    genotypes_asint = TargeneCore.get_genotype_encoding(param_files, ["RSID_2", "RSID_198"])
    @test genotypes_asint == false
    # Only Int is fine
    param_files = TargeneCore.load_param_files(joinpath("config", "param_1"))
    genotypes_asint = TargeneCore.get_genotype_encoding(param_files, ["RSID_2"])
    @test genotypes_asint == true
end

@testset "Test adjust_treatment_encoding!" begin
    ## With string encoded variants
    genotypes = DataFrame(
        SAMPLE_ID = [1, 2, 3],
        RSID_198 = ["GA", "GG", "AA"],
        RSID_2 = ["AG", "GG", "AA"],
    )
    # AG is not in the genotypes bu GA is
    param_file = TargeneCore.load_param_files(joinpath("config", "param_2_with"))[1]
    @test param_file["Parameters"][1]["RSID_198"] == Dict("case" => "AG", "control" => "AA")
    @test param_file["Parameters"][2]["RSID_198"] == "AG"
    TargeneCore.adjust_treatment_encoding!(param_file, genotypes)
    @test param_file["Parameters"][1]["RSID_198"] == Dict("case" => "GA", "control" => "AA")
    @test param_file["Parameters"][2]["RSID_198"] == "GA"
    # If the allele is not present 
    param_file = TargeneCore.load_param_files(joinpath("config", "param_2_with"))[1]
    genotypes.RSID_198 = ["AA", "GG", "AA"]
    @test_throws TargeneCore.AbsentAlleleError("RSID_198", "AG") TargeneCore.adjust_treatment_encoding!(param_file, genotypes)
    
    ## With integer encoded variants
    genotypes = DataFrame(
        SAMPLE_ID = [1, 2, 3],
        RSID_2 = [1, 2, 0],
    )
    param_file = TargeneCore.load_param_files(joinpath("config", "param_1.yaml"))[1]
    copied_param_file = deepcopy(param_file)
    TargeneCore.adjust_treatment_encoding!(param_file, genotypes)
    @test param_file == copied_param_file
    genotypes.RSID_2 = [2, 2, 0]
    @test_throws TargeneCore.AbsentAlleleError("RSID_2", 1) TargeneCore.adjust_treatment_encoding!(param_file, genotypes)

end
#####################################################################
###############           END-TO-END TESTS            ###############
#####################################################################


@testset "Test tmle_inputs from-param-files: scenario 1" begin
    # Genotypes encoded as strings
    parsed_args = Dict(
        "from-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_2_with_string_case_control.yaml")), 
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

    # Data Files
    data = CSV.read("final.data.csv", DataFrame)
    @test names(data) == [
        "SAMPLE_ID", "BINARY_1", "BINARY_2", "CONTINUOUS_1", 
        "CONTINUOUS_2", "COV_1", "21003", "22001", "TREAT_1", 
        "PC1", "PC2", "RSID_2", "RSID_198"
    ]
    # Parameter file:
    input_param_file = YAML.load_file(parsed_args["from-param-files"]["param-prefix"])
    out_param_file = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 1))
    # Added sections, all but 22001 are considered traits
    @test out_param_file["W"] == ["PC1", "PC2"]
    @test out_param_file["C"] == [22001]
    @test out_param_file["Y"] == ["BINARY_1", "BINARY_2", "CONTINUOUS_1", "CONTINUOUS_2", "COV_1", "21003", "TREAT_1"]
    # The rest is unchanged
    @test out_param_file["T"] == input_param_file["T"]
    @test out_param_file["Parameters"] == input_param_file["Parameters"]

    cleanup()
end

@testset "Test tmle_inputs from-param-files: scenario 2" begin
    parsed_args = Dict(
        "from-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_2.yaml")), 
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
    iate_param = Dict[Dict("name"=>"IATE", "TREAT_1"=>Dict("case"=>1, "control"=>0), "RSID_2"=>Dict("case"=>1, "control"=>0))]
    # 2 first variables in a batch
    param_file_1 = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 1))
    @test param_file_1 == Dict(
        "Y"          => ["BINARY_1", "BINARY_2"],
        "W"          => ["PC1", "PC2"],
        "T"          => ["RSID_2", "TREAT_1"],
        "Parameters" => iate_param
    )
    # CONTINUOUS_1 and CONTINUOUS_2
    param_file_2 = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 2))
    @test param_file_2 == Dict(
        "Y"          => ["CONTINUOUS_1", "CONTINUOUS_2"],
        "W"          => ["PC1", "PC2"],
        "T"          => ["RSID_2", "TREAT_1"],
        "Parameters" => iate_param
    )
    # COV_1 and 21003
    param_file_3 = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 3))
    @test param_file_3 == Dict(
        "Y"          => ["COV_1", "21003"],
        "W"          => ["PC1", "PC2"],
        "T"          => ["RSID_2", "TREAT_1"],
        "Parameters" => iate_param
    )
    # 22001
    param_file_4 = YAML.load_file(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 4))
    @test param_file_4 == Dict(
        "Y"          => ["22001"],
        "W"          => ["PC1", "PC2"],
        "T"          => ["RSID_2", "TREAT_1"],
        "Parameters" => iate_param
    )
    @test !isfile(TargeneCore.yaml_out_path(parsed_args["out-prefix"], 5))
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
        "positivity-constraint" => 0.,
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