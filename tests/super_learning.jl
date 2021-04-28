using GenesInteraction
using MLJ, Test, DataFrames, CategoricalArrays, MLJLinearModels


X = DataFrame(
    locus_1 = categorical(['A', 'C', 'C', 'G', 'A', 'C']),
    locus_2 = categorical(['C', 'T', 'T', 'T', 'C', 'C']),
    age = [45., 56., 34., 65., 24., 54.]
)

y = categorical([0, 0, 0, 1, 1, 0])

library = [@pipeline(OneHotEncoder(), 
                    LogisticClassifier(), 
                    name="MyPipeline")]

GenesInteraction.superlearn(library, X, y;nfolds=2)