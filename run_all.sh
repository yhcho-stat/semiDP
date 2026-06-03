#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

PYTHON_BIN="${PYTHON:-python3}"

Rscript scripts/semidp_gaussian_uniform.R
Rscript scripts/semidp_gaussian_increasing.R
Rscript scripts/semidp_Knorm_uniform.R
Rscript scripts/semidp_Knorm_linear.R
"${PYTHON_BIN}" scripts/make_paper_plots.py
Rscript scripts/sim_grouped_regression_mu_gdp.R
