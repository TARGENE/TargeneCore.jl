function cleanup(;prefix="final.")
    for file in readdir()
        if startswith(file, prefix)
            rm(file)
        end
    end
end


function make_estimands_configuration()
    estimands = [
        IATE(
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
        CM(
            outcome = "ALL",
            treatment_values = (RSID_2 = "GG", RSID_198 = "GA"),
            treatment_confounders = (RSID_2 = [], RSID_198 = []),
            outcome_extra_covariates = [22001]
        )
    ]
    return Configuration(estimands=estimands)
end

function make_estimands_configuration_no_wildcard()
    estimands = [
    IATE(
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