abstract type TMLEstimator end

pvalue(tmle::TMLEstimator) = 2*(1 - cdf(Normal(0, 1), abs(tmle.estimate/tmle.std_error)))

confint(tmle::TMLEstimator) = (tmle.estimate - 1.96tmle.stderror, tmle.estimate + 1.96tmle.stderror)