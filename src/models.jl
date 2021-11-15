###############################################################################
# REGRESSORS

SKLinearRegressor = @load LinearRegressor pkg=ScikitLearn verbosity=0
KNNRegressor = @load KNNRegressor pkg=NearestNeighborModels verbosity=0
LinearRegressor = @load LinearRegressor pkg=MLJLinearModels verbosity=0
XGBoostRegressor = @load XGBoostRegressor pkg=XGBoost verbosity=0

###############################################################################
# CLASSIFIERS

SKLogisticClassifier = @load LogisticClassifier pkg=ScikitLearn verbosity=0
LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels verbosity=0
KNNClassifier = @load KNNClassifier pkg=NearestNeighborModels verbosity=0
XGBoostClassifier = @load XGBoostClassifier pkg=XGBoost verbosity=0