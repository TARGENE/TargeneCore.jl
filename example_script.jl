using TargeneCore

testdir = joinpath(pkgdir(TargeneCore), "test")
cd(testdir)

include(joinpath(testdir, "runtests.jl"))