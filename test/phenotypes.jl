module TestPhenotypes

using Test
using TargeneCore
using CSV
using DataFrames

@testset "Test tmle_phenotypes_batches" begin
    n_phenotypes = 10
    temp_file = "temp_phenotypes.csv"
    phenotypes_names = ["x$i" for i in 1:n_phenotypes]
    CSV.write(
        temp_file, 
        DataFrame(rand(10, n_phenotypes+1), vcat("SAMPLE_ID", phenotypes_names))
        )

    # Test basic mode
    parsed_args = Dict(
        "phenotypes-file" => temp_file,
        "batch-size" => "3"
    )
    TargeneCore.tmle_phenotypes_batches(parsed_args)

    expected_batches = [
        ["x1", "x2", "x3"], 
        ["x4", "x5", "x6"], 
        ["x7", "x8", "x9"], 
        ["x10"]
    ]
    for batch in 1:4
        batch_file = TargeneCore.batch_filename(batch)
        @test readlines(open(batch_file)) == expected_batches[batch]
        rm(batch_file)
    end

    # Test max mode
    parsed_args = Dict(
        "phenotypes-file" => temp_file,
        "batch-size" => "max"
    )
    TargeneCore.tmle_phenotypes_batches(parsed_args)

    batch_file = TargeneCore.batch_filename(1)
    @test readlines(open(batch_file)) == phenotypes_names
    rm(batch_file)
    rm(temp_file)

end

end;

true