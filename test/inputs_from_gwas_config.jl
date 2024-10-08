module TestGwasEstimands

using Test
using SnpArrays
using TargeneCore
using Arrow
using DataFrames
using Serialization
using TMLE
using CSV

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "testutils.jl"))

function get_summary_stats(estimands)
    outcomes = [TargeneCore.get_outcome(Ψ) for Ψ in estimands]
    results = DataFrame(ESTIMAND = estimands, OUTCOME = outcomes)
    return sort(combine(groupby(results, :OUTCOME), nrow), :OUTCOME)
end

function check_estimands_levels_order(estimands, snp_info)
    for Ψ in estimands
        # If the two components are present, the first is the 0 -> 1 and the second is the 1 -> 2
        variant = only(keys(Ψ.args[1].treatment_values))
        variant_info = filter(:snpid=>x->x==String(variant),snp_info)
        allele1, allele2 = variant_info.allele1[1], variant_info.allele2[1]
        
        # Here, we check if the order is sufficient to be able to compute non-linear effects any of these combinations will do
        if length(Ψ.args) == 2
            @test (Ψ.args[1].treatment_values[variant] == (control = allele1*allele1, case = allele1*allele2) && 
            Ψ.args[2].treatment_values[variant] == (control = allele1*allele2, case = allele2*allele2)) || 
            (Ψ.args[1].treatment_values[variant] == (control = allele2*allele2, case = allele1*allele2) && 
            Ψ.args[2].treatment_values[variant] == (control = allele1*allele2, case = allele1*allele1))
        else
            # Otherwise we check they are one or the other
            arg = only(Ψ.args)
            @test arg.treatment_values[variant] == (control = allele1*allele1, case = allele1*allele2)  ||
            arg.treatment_values[variant] == (control = allele2*allele2, case = allele1*allele2) ||
            arg.treatment_values[variant] == (control = allele1*allele2, case = allele2*allele2) ||
            arg.treatment_values[variant] == (control = allele1*allele2, case = allele1*allele1) 

        end
   end
end
@testset "Test inputs_from_config gwas: no positivity constraint" begin
    tmpdir = mktempdir()
    copy!(ARGS, [
        "estimation-inputs",
        joinpath(TESTDIR, "data", "config_gwas.yaml"),
        string("--traits-file=", joinpath(TESTDIR, "data", "ukbb_traits.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "ukbb_pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "genotypes" , "ukbb_1.")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--batchsize=5",
        "--verbosity=0",
        "--positivity-constraint=0"
    ])
    TargeneCore.julia_main()

    # Define SNP information to check string allele defintions
    snpdata = read_bed_chromosome(joinpath(TESTDIR, "data", "ukbb", "genotypes" , "ukbb_1."))
    snp_info = select(DataFrame(snpdata.snp_info), [:snpid, :allele1, :allele2])
    # Check dataset
    dataset = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(dataset) == (1940, 886)

    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    @test all(e isa JointEstimand for e in estimands)

    # There are 875 variants in the dataset
    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:BINARY_1, :BINARY_2, :CONTINUOUS_1, :CONTINUOUS_2, :TREAT_1], 
        nrow = repeat([875], 5)
    )

    check_estimands_levels_order(estimands, snp_info)
end

@testset "Test inputs_from_config gwas: positivity constraint" begin
    tmpdir = mktempdir()
    copy!(ARGS, [
        "estimation-inputs",
        joinpath(TESTDIR, "data", "config_gwas.yaml"),
        string("--traits-file=", joinpath(TESTDIR, "data", "ukbb_traits.csv")),
        string("--pcs-file=", joinpath(TESTDIR, "data", "ukbb_pcs.csv")),
        string("--genotypes-prefix=", joinpath(TESTDIR, "data", "ukbb", "genotypes" , "ukbb_1.")),
        string("--outprefix=", joinpath(tmpdir, "final")), 
        "--batchsize=5",
        "--verbosity=0",
        "--positivity-constraint=0.2"
    ])
    TargeneCore.julia_main()
    # Check dataset
    dataset = DataFrame(Arrow.Table(joinpath(tmpdir, "final.data.arrow")))
    @test size(dataset) == (1940, 886)
    # Check estimands
    estimands = []
    for file in readdir(tmpdir, join=true)
        if endswith(file, "jls")
            append!(estimands, deserialize(file).estimands)
        end
    end
    # The positivity constraint reduces the number of variants
    @test all(e isa JointEstimand for e in estimands)
    summary_stats = get_summary_stats(estimands)
    @test summary_stats == DataFrame(
        OUTCOME = [:BINARY_1, :BINARY_2, :CONTINUOUS_1, :CONTINUOUS_2, :TREAT_1], 
        nrow = repeat([777], 5)
    )
   
    check_estimands_levels_order(estimands)
end

end

true