module TestDataset

using Test
using TargeneCore
using Arrow
using DataFrames

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

@testset "Test make_dataset" begin
    tmpdir = mktempdir()
    out = joinpath(tmpdir, "dataset.arrow")
    copy!(ARGS, [
        "make-dataset",
        string("--traits-file=", joinpath(TESTDIR, "data", "traits_1.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "imputed" ,"ukbb")),
        string("--variants-file=", joinpath(TESTDIR, "data", "variants_list.txt")),
        string("--out=", out),
        "--call-threshold=0.8",
    ])
    TargeneCore.julia_main()

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