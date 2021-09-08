###############################################################################
# REGRESSORS

RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree verbosity = 0
NuSVR = @load NuSVR pkg=LIBSVM verbosity = 0
LinearRegressor = @load LinearRegressor pkg=MLJLinearModels verbosity = 0
KNNRegressor = @load KNNRegressor pkg=NearestNeighborModels verbosity = 0
NeuralNetworkRegressor = @load NeuralNetworkRegressor pkg=MLJFlux verbosity = 0
EvoTreeRegressor = @load EvoTreeRegressor pkg=EvoTrees verbosity = 0

###############################################################################
# CLASSIFIERS

RandomForestClassifier= @load RandomForestClassifier pkg=DecisionTree verbosity = 0
NuSVC = @load NuSVC pkg=LIBSVM verbosity = 0
LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity = 0
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels verbosity = 0
NeuralNetworkClassifier = @load NeuralNetworkClassifier pkg=MLJFlux verbosity = 0
EvoTreeClassifier = @load EvoTreeClassifier pkg=EvoTrees verbosity = 0