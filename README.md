# Formal Privacy Guarantees with Invariant Statistics: Reproducibility Code

This repository contains the simulation code and generated artifacts for the
numerical illustrations in the JMLR revision of **"Formal Privacy Guarantees
with Invariant Statistics."**

Generated CSV files and PNG figures are intentionally included so that the
figures can be inspected without rerunning the simulations.

## Directory Layout

```text
.
├── scripts/
│   ├── semidp_gaussian_uniform.R
│   ├── semidp_gaussian_increasing.R
│   ├── semidp_Knorm_uniform.R
│   ├── semidp_Knorm_linear.R
│   ├── make_paper_plots.py
│   ├── sim_grouped_regression_mu_gdp.R
│   └── sim_grouped_anova_mu_gdp.R
├── output/
│   ├── gaussian_uniform.csv
│   ├── gaussian_linear.csv
│   ├── Knorm_uniform.csv
│   └── Knorm_linear.csv
├── results/
│   ├── grouped_regression_mu_gdp.csv
│   └── grouped_anova_mu_gdp.csv
├── plots/
│   ├── gaussian_comparison.png
│   ├── Knorm_comparison.png
│   ├── grouped_regression_mu_gdp.png
│   ├── grouped_anova_mu_gdp.png
│   └── fig_7by6table.png
└── run_all.sh
```

## Why Both R and Python?

The code reflects the workflow used during the paper revision.

- The original contingency-table simulations were written in R. Those scripts
  generate CSV files in `output/`.
- The original contingency-table plots were assembled from those CSV files in a
  Python/Jupyter notebook. For a cleaner command-line workflow, the plotting
  code has been moved into `scripts/make_paper_plots.py`, which writes
  `plots/gaussian_comparison.png` and `plots/Knorm_comparison.png`.
- The newer grouped regression and group-means simulations were written as
  standalone R scripts. Each of these scripts writes both its CSV file and its
  PNG figure.

## Requirements

R packages:

- `MASS`
- `pracma`
- `lpSolve`
- `LaplacesDemon`

Python packages:

- `pandas`
- `matplotlib`

On macOS, the system `python3` may not include these packages. If you use
Anaconda or another environment, call that Python explicitly, for example:

```bash
/opt/anaconda3/bin/python scripts/make_paper_plots.py
```

The grouped regression and grouped ANOVA scripts use base R only.

## Reproducing the Figures

Run commands from the repository root.

### Contingency-table Gaussian simulations

```bash
Rscript scripts/semidp_gaussian_uniform.R
Rscript scripts/semidp_gaussian_increasing.R
```

Outputs:

- `output/gaussian_uniform.csv`
- `output/gaussian_linear.csv`

### Contingency-table K-norm simulations

```bash
Rscript scripts/semidp_Knorm_uniform.R
Rscript scripts/semidp_Knorm_linear.R
```

Outputs:

- `output/Knorm_uniform.csv`
- `output/Knorm_linear.csv`

The K-norm scripts use a rejection-sampling step and can take longer than the
Gaussian scripts.

### Legacy contingency-table plots

After the contingency CSV files exist, run:

```bash
python3 scripts/make_paper_plots.py
```

Outputs:

- `plots/gaussian_comparison.png`
- `plots/Knorm_comparison.png`

### Grouped regression with count invariants

```bash
Rscript scripts/sim_grouped_regression_mu_gdp.R
```

Outputs:

- `results/grouped_regression_mu_gdp.csv`
- `plots/grouped_regression_mu_gdp.png`

### Group means / one-way ANOVA with count invariants

```bash
Rscript scripts/sim_grouped_anova_mu_gdp.R
```

Outputs:

- `results/grouped_anova_mu_gdp.csv`
- `plots/grouped_anova_mu_gdp.png`

## One-command Reproduction

To rerun all simulations and plotting scripts:

```bash
bash run_all.sh
```

This reruns the K-norm simulations as well, so it may take a while.
If your plotting packages are installed in a non-default Python, set `PYTHON`:

```bash
PYTHON=/opt/anaconda3/bin/python bash run_all.sh
```

## Notes on Included Artifacts

- `plots/fig_7by6table.png` is the illustrative sensitivity-space table figure
  used in the manuscript.
- `output/` stores legacy contingency-table simulation summaries.
- `results/` stores the newer grouped regression and group-means simulation
  summaries.
- `plots/` stores all generated PNG figures used by the manuscript.
