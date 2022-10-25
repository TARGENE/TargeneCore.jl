module TestConfoundersPreparation

using Test
using SnpArrays
using TargeneCore
using DataFrames
using CSV
using BGEN

function clean(file)
    if endswith(file, ".bgen")
        out_samples = file[1:end-4]*"sample"
        rm(file)
        rm(out_samples)
    elseif endswith(file, ".bed")
        rm(file)
        rm(file[1:end-3]*"bim")
        rm(file[1:end-3]*"fam")
    end
end

@testset "Various functions" begin
    # Test issnp
    @test TargeneCore.issnp("A") == true
    @test TargeneCore.issnp("AA") == false
    
    # Test bounds
    @test TargeneCore.bounds(0, -10, 10) == (-10^4, 10^4)
    @test TargeneCore.bounds(15000, 1000, 30000) == (1000, 30000)

    # Test ld_blocks_by_chr
    lb_blocks_by_chr = TargeneCore.ld_blocks_by_chr(joinpath("data", "VDR_LD_blocks.txt"))
    @test lb_blocks_by_chr[1][1, :lower_bound] == -9990.0
    @test lb_blocks_by_chr[1][1, :upper_bound] == 10010.0
    @test lb_blocks_by_chr[1][1, :chr] == 3

    @test lb_blocks_by_chr[2][1, :lower_bound] == 4.8262895e7
    @test lb_blocks_by_chr[2][1, :upper_bound] == 4.8282895e7
    @test lb_blocks_by_chr[2][1, :chr] == 12

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

@testset "Test filter_genetic_file: bgen & bed" begin
    # Scenario 1: same results BGEN / BED
    parsed_args = Dict(
        "input"  => joinpath("data", "bgen", "example.8bits.bgen"),
        "output" => joinpath("data", "filtered-example.8bits.bgen"),
        "qcfile" => joinpath("data", "qc_file_example.8bits.csv"),
        "ld-blocks" => joinpath("data", "LD_blocks.example.8bits.txt"),
        "maf-threshold" => 0.22,
        "traits" => joinpath("data", "sample_ids.example.8bits.txt")
    )
    # RSID_101: not in both arrays
    # RSID_2: not in both arrays
    # RSID_102: not all batches pass qc
    # RSID_3: pass
    # RSID_103: not a SNP
    # RSID_4: maf under threshold
    # RSID_104: no in both arrays
    # RSID_5: pass
    # RSID_105: in LD
    expected_snps = ["RSID_3", "RSID_5"]
    expected_person = ["sample_001", "sample_002", "sample_003", "sample_004", "sample_005"]
    filter_ukb_genetic_file(parsed_args)
    out_samples = parsed_args["output"][1:end-4]*"sample"
    out_bgen = Bgen(parsed_args["output"]; sample_path=out_samples)
    @test [x.rsid for x in BGEN.iterator(out_bgen)] == expected_snps
    out_bgen.samples == expected_person
    clean(parsed_args["output"])

    parsed_args["input"] = joinpath("data", "bed", "example.8bits.bed")
    parsed_args["output"] = joinpath("data", "filtered-example.8bits.bed")
    filter_ukb_genetic_file(parsed_args)
    snp_data = SnpData(parsed_args["output"][1:end-4])
    @test snp_data.snp_info.snpid == expected_snps
    @test snp_data.person_info.iid == expected_person
    clean(parsed_args["output"])

    # Scenario 2: same results BGEN / BED
    # Now RSID_4 will pass
    parsed_args["maf-threshold"] = 0.01
    expected_snps = ["RSID_3", "RSID_4", "RSID_5"]
    filter_ukb_genetic_file(parsed_args)
    snp_data = SnpData(parsed_args["output"][1:end-4])
    @test snp_data.snp_info.snpid == expected_snps
    @test snp_data.person_info.iid == expected_person
    clean(parsed_args["output"])

    parsed_args["input"] = joinpath("data", "bgen", "example.8bits.bgen")
    parsed_args["output"] = joinpath("data", "filtered-example.8bits.bgen")

    filter_ukb_genetic_file(parsed_args)
    out_samples = parsed_args["output"][1:end-4]*"sample"
    out_bgen = Bgen(parsed_args["output"]; sample_path=out_samples)
    @test [x.rsid for x in BGEN.iterator(out_bgen)] == expected_snps
    out_bgen.samples == expected_person
    clean(parsed_args["output"])

end


@testset "Test merge_beds" begin
    genotypes_dir = joinpath("data", "ukbb", "genotypes")
    chr_prefix = joinpath(genotypes_dir, "ukbb")

    parsed_args = Dict(
            "input"  => chr_prefix,
            "output" => joinpath(genotypes_dir, "ukbb_merged")
        )
    merge_beds(parsed_args)
    merged = SnpData(parsed_args["output"])

    @test length(unique(merged.snp_info.chromosome)) == 3
    # Clean

    for ext in [".bed", ".bim", ".fam"]
        rm(parsed_args["output"]*ext)
    end

end


@testset "Test adapt_flashpca" begin
    parsed_args = Dict(
        "input" => joinpath("data", "flashpca_output.txt"),
        "output" => "flashpca_after_adapt.csv",
    )
    adapt_flashpca(parsed_args)

    output = CSV.File(parsed_args["output"]) |> DataFrame

    @test names(output) == ["SAMPLE_ID",
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"]
    @test size(output) == (2, 7)

    rm(parsed_args["output"])
end

end;

true