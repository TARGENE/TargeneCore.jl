module TestDataset

using Test
using TargeneCore
using Arrow
using DataFrames

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

@testset "Test generate_dataset(parsed_args)" begin
    tmpdir = mktempdir()
    out = joinpath(tmpdir, "dataset.arrow")
    parsed_args = Dict(
        "verbosity" => 0,
        "traits" => joinpath(TESTDIR, "data", "traits_1.csv"),
        "confounders" => joinpath(TESTDIR, "data", "pcs.csv"),
        "bgen-prefix" => joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb"),
        "variants" => joinpath(TESTDIR, "data", "variants_list.txt"),
        "call-threshold" => 0.8,
        "out" => out
    )
    generate_dataset(parsed_args)
    dataset = Arrow.Table(out) |> DataFrame
    @test names(dataset) == [
        "SAMPLE_ID",
        "BINARY_1",
        "BINARY_2",
        "CONTINUOUS_1",
        "CONTINUOUS_2",
        "COV_1",
        "21003",
        "22001",
        "TREAT_1",
        "PC1",
        "PC2",
        "RSID_103",
        "RSID_198"
    ]
    @test size(dataset, 1) == 490

end

end

true