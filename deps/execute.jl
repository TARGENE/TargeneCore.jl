using TargeneCore

@info "Running precompilation script."
# Run workload
TEST_DIR = joinpath(pkgdir(TargeneCore), "test")
push!(LOAD_PATH, TEST_DIR)
include(joinpath(TEST_DIR, "runtests.jl"))