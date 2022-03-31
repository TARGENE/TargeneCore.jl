using TMLEEpistasis

testdir = joinpath(pkgdir(TMLEEpistasis), "test")
cd(testdir)
include(joinpath(testdir, "runtests.jl"))