using TargeneCore

TESTDIR = joinpath(pkgdir(TargeneCore), "test")


include(joinpath(TESTDIR, "confounders.jl"))
include(joinpath(TESTDIR, "tmle_inputs", "tmle_inputs.jl"))
include(joinpath(TESTDIR,"tmle_inputs", "from_actors.jl"))
include(joinpath(TESTDIR, "tmle_inputs", "from_param_files.jl"))