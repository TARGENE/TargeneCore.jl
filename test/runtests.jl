using TargeneCore
using Test

TESTDIR = joinpath(pkgdir(TargeneCore), "test")

@test include(joinpath(TESTDIR, "utils.jl"))
@test include(joinpath(TESTDIR, "confounders.jl"))
@test include(joinpath(TESTDIR, "tl_inputs", "tl_inputs.jl"))
@test include(joinpath(TESTDIR, "tl_inputs", "from_actors.jl"))
@test include(joinpath(TESTDIR, "tl_inputs", "from_param_files.jl"))
@test include(joinpath(TESTDIR, "tl_inputs", "allele_independent_estimands.jl"))
@test include(joinpath(TESTDIR, "tl_inputs", "permutation_test.jl"))
@test include(joinpath(TESTDIR, "tl_inputs", "random_variants_test.jl"))