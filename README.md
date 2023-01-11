# TargeneCore.jl

![GitHub Workflow Status (branch)](https://img.shields.io/github/workflow/status/TARGENE/TargeneCore.jl/CI/main?label=Build%20main)
![Codecov branch](https://img.shields.io/codecov/c/github/TARGENE/TargeneCore.jl/main?label=Coverage%20main)
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/TARGENE/TargeneCore.jl)

This project provides various functionalities that are used thoughout the [targene pipeline](https://github.com/TARGENE/targene-pipeline).

## Preparation of TMLE input data

The associated script can be used by:

```bash
julia --project --startup-file=no tmle_inputs.jl --help
```

The purpose of this step is to generate a set of data and [parameters configuration files](test/config) that can be further used by the [Targeted Estimation executable](https://github.com/TARGENE/TargetedEstimation.jl).

The data files generation process works as follows. All input data files are merged together and further divided as:

- Genetic & Extra confounders &rarr; Confounders
- SNPs & Extra treatments &rarr; Treatments
- Binary phenotypes &rarr; Binary phenotypes
- Continuous phenotypes &rarr; Continuous phenotypes
- Covariates &rarr; Covariates

Apart from the SNP data, all input data files are supposed to be given as CSV filles with a `SAMPLE_ID` column that will be used for merging. There are two ways by which parameters configuration files and SNP data are generated.

### Provided configuration files

For this mode, use the `with-param-files` command and `--param-prefix` option. The SNPs of interest are read from those configuration files that will be validated against the actual data to ensure correctness.

### MechanismActors strategy

For this strategy, use the `with-actors` command. and `--asb-prefix` and `--trans-actors` options. This will generate pairwise interaction parameters between bQTLs output by the [Baal-ChIP pipeline](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf) and trans-actors eQTLs from a CSV file. Additionally, if template parameters configuration files containing extra treatments are provided, nth-order interaction parameters will be generated.
