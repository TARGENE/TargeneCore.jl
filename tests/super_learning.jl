using GenesInteraction
using Distributions, MLJ, Test, DataFrames, CategoricalArrays, MLJLinearModels

function simulation_dataset(;n=10000)
    Age = rand(DiscreteUniform(30, 70), n)
    X_1 = 1 ./ (1 .+ exp.(-0.005*Age .- 0.01)) .< rand(Uniform(0,1), n)
    X_2 = 1 ./ (1 .+ exp.(-0.01*Age)) .< rand(Uniform(0,1), n)
    X = DataFrame(
        locus_1 = categorical(X_1),
        locus_2 = categorical(X_2),
        age = Age
    )

    y = 1 ./ (1 .+ exp.(-0.01*Age - 0.94*X_1 + 0.4*X_2)) .< rand(Uniform(0,1), n)
    X, categorical(y)
end


function simple_library()
    library = []
    for λ in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1., 10., 100., 1000.]
        lr = LogisticClassifier(lambda=λ)
        push!(library, lr)
    end
    return library
end

lib = simple_library()

X, y = simulation_dataset(;n=1000000)

sl = SuperLearner(lib, LogisticClassifier(), 3, false)

pipe = @pipeline(OneHotEncoder(;drop_last=true), sl)

mach = machine(pipe, X, y)

fit!(mach)

# I'd like to access the different machines in a user friendly manner to test the implementation
# Not sure, but it seems that using the "report" attribute is the way to go
sl_machine = mach.report.machines[2]

# This is a 37 list of machines
sl_sub_machines = sl_machine.report.machines

# For instance, what is this machine corresponding to?
# I'd like to know if this is the global fit or a i-th fold fit
# I guess the order of creation is maintained in this list so I could probably
# write a method to sort this out but a naming mechanism would be more elegant in my view
sub_machine = sl_sub_machines[10]

# For completeness, I will probably be interested in optionally evaluating folds for ech machine 
# and being able to retrieve machines scores on each fold, and the average score.


@pipeline(OneHotEncoder(;drop_last=true), sl)