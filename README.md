# TMLEEpistasis.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://olivierlabayle.github.io/TMLEEpistasis.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://olivierlabayle.github.io/TMLEEpistasis.jl/dev)
[![Build Status](https://github.com/olivierlabayle/TMLEEpistasis.jl/workflows/CI/badge.svg)](https://github.com/olivierlabayle/TMLEEpistasis.jl/actions)
[![Coverage](https://codecov.io/gh/olivierlabayle/TMLEEpistasis.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/olivierlabayle/TMLEEpistasis.jl)


The purpose of this project is to provide a mean for the estimation of effect sizes of potentially interacting variants using the Targeted Learning framework.

## Prerequisites

This project has important non Julia dependencies that may be difficult to install. The easiest way to use the facilities offered here is to use the dedicated [docker image](https://hub.docker.com/repository/docker/olivierlabayle/tmle-epistasis).

## Description of the command line interface

From the project's root directory:

```bash
julia --project ukbb.jl --help
```
