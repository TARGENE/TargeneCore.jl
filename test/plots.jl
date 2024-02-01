module TestPlots

using TargeneCore
using Test
using CairoMakie

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "testutils.jl"))

@testset "Test generate_summary_plots" begin
    tmpdir = mktempdir()
    prefix = joinpath(tmpdir, "tmle_output")
    make_fake_outputs(; prefix=prefix)
    parsed_args = Dict(
        "verbosity" => 0,
        "results-prefix" => string(prefix, ".hdf5"),
        "out-prefix" => prefix
    )
    fig = generate_summary_plots(parsed_args)
    @test length(fig.content) == 2
    @test fig.content[1] isa CairoMakie.Axis
    @test fig.content[2] isa CairoMakie.Legend
    @test isfile(TargeneCore.filepath_from_prefix(prefix; filename="QQ.png"))
end

end

true