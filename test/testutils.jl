using TMLE
using TMLECLI

function make_dataset(;n=100, rng=StableRNG(123))
    return DataFrame(
        "rs10043934" => rand(rng, ["GA", "GG", "AA"], n),
        "rs117913124" => rand(rng, ["GA", "GG", "AA"], n),
        "RSID_103" => rand(rng, ["GA", "GG", "AA"], n),
        "RSID_104" => rand(rng, ["GA", "GG", "AA"], n),
        "High light scatter reticulocyte percentage" => ones(n),
        "L50-L54 Urticaria and erythema" =>  ones(n)
    )
end

function make_estimands()
    AIE₁ = AIE(
        outcome = "High light scatter reticulocyte percentage",
        treatment_values = (
            rs10043934 = (case="GA", control="GG"), 
            RSID_103 = (case="GA", control="GG")
        ),
        treatment_confounders = (
            rs10043934 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6), 
            RSID_103 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6)
        ),
        outcome_extra_covariates = ("Age-Assessment", "Genetic-Sex")
    )
    AIE₂ = AIE(
        outcome = "High light scatter reticulocyte percentage",
        treatment_values = (
            rs10043934 = (case="GA", control="GG"), 
            RSID_103 = (case="AA", control="GA")
        ),
        treatment_confounders = (
            rs10043934 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6), 
            RSID_103 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6)
        ),
        outcome_extra_covariates = ("Age-Assessment", "Genetic-Sex")
    )
    ATE₁ = ATE(
        outcome = "L50-L54 Urticaria and erythema",
        treatment_values = (
            rs117913124 = (case="GA", control="GG"), 
            RSID_104 = (case="GA", control="GG")
        ),
        treatment_confounders = (
            rs117913124 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6), 
            RSID_104 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6)
        ),
        outcome_extra_covariates = ("Age-Assessment", "Genetic-Sex")
    )
    jointAIE = JointEstimand(AIE₁, AIE₂)
    return (AIE₁, AIE₂, jointAIE, ATE₁)
end

function make_estimates()
    AIE₁, AIE₂, jointAIE, ATE₁ = make_estimands()
    AIE₁ = TMLE.TMLEstimate(
        estimand = AIE₁,
        estimate = -1.,
        std = 0.003,
        n = 10,
        IC = []
    )
    AIE₂ = TMLE.TMLEstimate(
        estimand = AIE₂,
        estimate = -0.003,
        std = 0.003,
        n = 10,
        IC = []
    )

    jointAIE = TMLE.JointEstimate(
        estimand = jointAIE,
        estimates = (AIE₁, AIE₂),
        cov = [
        0.003 0.
        0. 0.003
        ],
        n = 10
    )

    ATE₁ = TMLE.TMLEstimate(
        estimand = ATE₁,
        estimate = 0.003,
        std = 0.003,
        n = 20,
        IC = []
    )

    failed_estimate = TMLECLI.FailedEstimate(
        ATE(
            outcome = "L50-L54 Urticaria and erythema",
            treatment_values = (
                rs117913124 = (case="GA", control="GG"), 
                RSID_104 = (case="GA", control="GG")
            ),
            treatment_confounders = (
                rs117913124 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6), 
                RSID_104 = (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6)
            ),
            outcome_extra_covariates = ("Age-Assessment", "Genetic-Sex")
        ),
        "Could not fluctuate"
    )
    return [
        (TMLE=AIE₁, OSE=AIE₁), 
        (TMLE=AIE₂, OSE=AIE₂), 
        (TMLE=jointAIE, OSE=jointAIE), 
        (TMLE=ATE₁, OSE=ATE₁), 
        (TMLE=failed_estimate, OSE=failed_estimate)
    ]
end

function save(estimates; prefix="tmle_output")
    # Batch 1
    outputs = TMLECLI.Outputs(hdf5=prefix*"_1.hdf5")
    TMLECLI.initialize(outputs)
    TMLECLI.update(outputs, estimates[1:3])
    TMLECLI.finalize(outputs)
    # Batch 2
    outputs = TMLECLI.Outputs(hdf5=prefix*"_2.hdf5")
    TMLECLI.initialize(outputs)
    TMLECLI.update(outputs, estimates[4:end])
    TMLECLI.finalize(outputs)
end

make_fake_outputs(estimates_generator=make_estimates; prefix="tmle_output") = 
    save(estimates_generator(); prefix=prefix)

### Fixtures for inputs_from_estimands

function make_estimands_configuration()
    estimands = [
        AIE(
            outcome = "ALL",
            treatment_values = (RSID_2 = (case = "AA", control = "GG"), TREAT_1 = (case = 1, control = 0)),
            treatment_confounders = (RSID_2 = [], TREAT_1 = [])
        ),
        ATE(
            outcome = "ALL",
            treatment_values = (RSID_2 = (case = "AA", control = "GG"),),
            treatment_confounders = (RSID_2 = [22001], ),
            outcome_extra_covariates = ["COV_1", 21003]
        ),
        CM(
            outcome = "ALL",
            treatment_values = (RSID_2 = "AA", ),
            treatment_confounders = (RSID_2 = [22001],),
            outcome_extra_covariates = ["COV_1", 21003]
        ),
        ATE(
            outcome = "ALL",
            treatment_values = (RSID_2 = (case = "AA", control = "GG"), RSID_198 = (case = "AG", control = "AA")),
            treatment_confounders = (RSID_2 = [], RSID_198 = []),
            outcome_extra_covariates = [22001]
        ),
        JointEstimand(
            CM(
                outcome = "ALL",
                treatment_values = (RSID_2 = "GG", RSID_198 = "AG"),
                treatment_confounders = (RSID_2 = [:PC1], RSID_198 = [:PC2]),
                outcome_extra_covariates = [22001]
            ),
            CM(
                outcome = "ALL",
                treatment_values = (RSID_2 = "AA", RSID_198 = "AG"),
                treatment_confounders = (RSID_2 = [:PC1], RSID_198 = [:PC2]),
                outcome_extra_covariates = [22001]
            )
        )
    ]
    return Configuration(estimands=estimands)
end

function make_estimands_configuration_no_wildcard()
    estimands = [
    AIE(
        outcome = "BINARY_1",
        treatment_values = (RSID_2 = (case = "AA", control = "GG"), TREAT_1 = (case = 1, control = 0)),
        treatment_confounders = (RSID_2 = [], TREAT_1 = [])
    ),
    ATE(
        outcome =  "CONTINUOUS_2",
        treatment_values = (RSID_2 = (case = "AA", control = "GG"),),
        treatment_confounders = (RSID_2 = [22001], ),
        outcome_extra_covariates = ["COV_1", 21003]
    ),
    CM(
        outcome =  "ALL",
        treatment_values = (RSID_2 = "AA", ),
        treatment_confounders = (RSID_2 = [22001],),
        outcome_extra_covariates = ["COV_1", 21003]
    )
    ]
    return Configuration(estimands=estimands)
end

function make_estimands_configuration_file(config_generator=make_estimands_configuration)
    dir = mktempdir()
    config = config_generator()
    filename = joinpath(dir, "configuration.yaml")
    TMLE.write_yaml(filename, config)
    return filename
end