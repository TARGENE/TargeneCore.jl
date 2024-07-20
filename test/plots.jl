module TestPlots

using TargeneCore
using Test
using CairoMakie

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

include(joinpath(TESTDIR, "testutils.jl"))

@testset "Test summary_plots" begin
    tmpdir = mktempdir()
    prefix = joinpath(tmpdir, "tmle_output")
    make_fake_outputs(; prefix=prefix)
    copy!(ARGS, [
        "summary-plots",
        string(prefix, ".hdf5"),
        string("--outprefix=", prefix)
    ])
    TargeneCore.julia_main()
    @test isfile(TargeneCore.filepath_from_prefix(prefix; filename="QQ.png"))
end

end

true