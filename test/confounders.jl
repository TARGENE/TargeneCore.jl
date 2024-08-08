module TestConfoundersPreparation

using Test
using SnpArrays
using TargeneCore
using DataFrames
using CSV

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

@testset "Various functions" begin
    # Test issnp
    @test TargeneCore.issnp("A") == true
    @test TargeneCore.issnp("AA") == false
    
    #Test bounds
    @test TargeneCore.bounds(0, -10, 10) == (-10^4, 10^4)
    @test TargeneCore.bounds(15000, 1000, 30000) == (1000, 30000)

    # Test notin_ldblocks
    ldblocks = DataFrame(chr        = [1, 3],
                        lower_bound = [10, 19],
                        upper_bound = [100, 109])
    ldblocks = groupby(ldblocks, :chr)

    snp_info = DataFrame(chromosome=[1, 1, 2, 3, 3],
                        position=[21, 5, 9, 102, 4])

    expected_snps = DataFrame(chromosome = [1, 2, 3],
                              position   = [5, 9, 4])
    @test filter(x -> TargeneCore.notin_ldblocks(x, ldblocks), snp_info) == expected_snps
end

@testset "Test filter_chromosome" begin
    tmpdir = mktempdir()
    output = joinpath(tmpdir, "filtered-mouse")
    # All options provided
    copy!(ARGS, [
        "filter-chromosome",
        SnpArrays.datadir("mouse"),
        output,
        joinpath(TESTDIR, "data", "sample_ids.txt"),
        string("--qc-file=", joinpath(TESTDIR, "data", "ukbb", "qcfile.txt")),
        string("--ld-blocks-file=", joinpath(TESTDIR, "data", "LD_blocks.txt")),
        "--maf-threshold=0.31"
    ])
    TargeneCore.julia_main()

    filtered = SnpData(output)

    @test filtered.snp_info.snpid == ["rs13476318"]
    @test size(filtered.snparray) == (5, 1)
    @test filtered.person_info.iid == 
        ["A048005080", "A048006063", "A048006555", "A048007096", "A048010273"]
    
    # No QC file provided
    copy!(ARGS, [
        "filter-chromosome",
        SnpArrays.datadir("mouse"),
        output,
        joinpath(TESTDIR, "data", "sample_ids.txt"),
        string("--ld-blocks-file=", joinpath(TESTDIR, "data", "LD_blocks.txt")),
        "--maf-threshold=0.495"
    ])
    TargeneCore.julia_main()

    filtered = SnpData(output)
    @test size(filtered.snparray) == (5, 88)
    @test filtered.person_info.iid == 
        ["A048005080", "A048006063", "A048006555", "A048007096", "A048010273"]
    
    # No ld-block file provided
    copy!(ARGS, [
        "filter-chromosome",
        SnpArrays.datadir("mouse"),
        output,
        joinpath(TESTDIR, "data", "sample_ids.txt"),
        "--maf-threshold=0.495"
    ])
    TargeneCore.julia_main()

    filtered = SnpData(output)
    ## More variants than the previous settings
    @test size(filtered.snparray) == (5, 95)
    @test filtered.person_info.iid == 
        ["A048005080", "A048006063", "A048006555", "A048007096", "A048010273"]
end

@testset "Test merge_beds" begin
    # Only works with relative paths at the moment see https://github.com/OpenMendel/SnpArrays.jl/issues/135
    genotypes_dir = relpath(joinpath(TESTDIR, "data", "ukbb", "genotypes"))
    chr_prefix = joinpath(genotypes_dir, "ukbb")
    tmpdir = mktempdir()
    output = joinpath(tmpdir, "ukbb_merged")
    copy!(ARGS, [
        "merge-beds",
        chr_prefix,
        output
    ])
    TargeneCore.julia_main()

    merged = SnpData(output)

    @test length(unique(merged.snp_info.chromosome)) == 3
end

@testset "Test adapt_flashpca" begin
    tmpdir = mktempdir()
    output_file = joinpath(tmpdir, "flashpca_after_adapt.csv")
    copy!(ARGS, [
        "adapt-flashpca",
        joinpath(TESTDIR, "data", "flashpca_output.txt"),
        output_file
    ])
    TargeneCore.julia_main()

    output = CSV.File(output_file) |> DataFrame

    @test names(output) == ["SAMPLE_ID",
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"]
    @test size(output) == (2, 7)

    rm(output_file)
end


end;

true
