using TMLEEpistasis

testdir = joinpath(pkgdir(TMLEEpistasis), "test")
cd(testdir)
include("models.jl")
include("estimators.jl")
include("utils.jl")
include(joinpath(testdir, "ukbb.jl"))