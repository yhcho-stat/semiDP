#!/usr/bin/env python3
"""Create manuscript plots from saved CSV files.

The legacy contingency-table simulations save CSV files from R. This script is
the script version of the original plotting notebook and writes the corresponding
PNG files to the plots/ directory.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
PLOTS.mkdir(parents=True, exist_ok=True)


def make_gaussian_plot() -> None:
    gaussian_linear_df = pd.read_csv(RESULTS / "gaussian_linear.csv")
    gaussian_uniform_df = pd.read_csv(RESULTS / "gaussian_uniform.csv")

    k_values = gaussian_linear_df["k"].tolist()
    loss_semi_dp_linear = gaussian_linear_df["LossSemiDP"].tolist()
    loss_naive_linear = gaussian_linear_df["LossNaive"].tolist()
    loss_semi_dp_uniform = gaussian_uniform_df["LossSemiDP"].tolist()
    loss_naive_uniform = gaussian_uniform_df["LossNaive"].tolist()

    plt.figure(figsize=(14, 6))

    plt.subplot(1, 2, 1)
    plt.plot(k_values, loss_semi_dp_linear, marker="o", linestyle="-", label="Our Mechanism")
    plt.plot(k_values, loss_naive_linear, marker="o", linestyle="--", label="Naive Mechanism")
    plt.xlabel("Size k")
    plt.ylabel("L2 Cost")
    plt.title("L2-Costs Comparison: Model I")
    plt.legend()
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(k_values, loss_semi_dp_uniform, marker="o", linestyle="-", label="Our Mechanism")
    plt.plot(k_values, loss_naive_uniform, marker="o", linestyle="--", label="Naive Mechanism")
    plt.xlabel("Size k")
    plt.ylabel("L2 Cost")
    plt.title("L2-Costs Comparison: Model II")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(PLOTS / "gaussian_comparison.png")
    plt.close()


def get_loss_data(df: pd.DataFrame, k_value: int, loss_column: str, epsilon_values: list[float]) -> list[float]:
    filtered_df = df[df["k"] == k_value]
    return [
        filtered_df.loc[filtered_df["epsilon"] == eps, loss_column].values[0]
        for eps in epsilon_values
    ]


def make_knorm_plot() -> None:
    knorm_linear_df = pd.read_csv(RESULTS / "Knorm_linear.csv")
    knorm_uniform_df = pd.read_csv(RESULTS / "Knorm_uniform.csv")

    epsilon_values = [0.1, 0.5, 1]

    model_1_l1 = [
        get_loss_data(knorm_linear_df, 2, "LossL1", epsilon_values),
        get_loss_data(knorm_linear_df, 3, "LossL1", epsilon_values),
    ]
    model_1_l2 = [
        get_loss_data(knorm_linear_df, 2, "LossL2", epsilon_values),
        get_loss_data(knorm_linear_df, 3, "LossL2", epsilon_values),
    ]
    model_1_linf = [
        get_loss_data(knorm_linear_df, 2, "LossLinf", epsilon_values),
        get_loss_data(knorm_linear_df, 3, "LossLinf", epsilon_values),
    ]
    model_1_knorm = [
        get_loss_data(knorm_linear_df, 2, "LossKnorm", epsilon_values),
        get_loss_data(knorm_linear_df, 3, "LossKnorm", epsilon_values),
    ]

    model_2_l1 = [
        get_loss_data(knorm_uniform_df, 2, "LossL1", epsilon_values),
        get_loss_data(knorm_uniform_df, 3, "LossL1", epsilon_values),
    ]
    model_2_l2 = [
        get_loss_data(knorm_uniform_df, 2, "LossL2", epsilon_values),
        get_loss_data(knorm_uniform_df, 3, "LossL2", epsilon_values),
    ]
    model_2_linf = [
        get_loss_data(knorm_uniform_df, 2, "LossLinf", epsilon_values),
        get_loss_data(knorm_uniform_df, 3, "LossLinf", epsilon_values),
    ]
    model_2_knorm = [
        get_loss_data(knorm_uniform_df, 2, "LossKnorm", epsilon_values),
        get_loss_data(knorm_uniform_df, 3, "LossKnorm", epsilon_values),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    titles = ["Model 1: k = 2", "Model 1: k = 3", "Model 2: k = 2", "Model 2: k = 3"]
    data_sets = [
        (model_1_l1[0], model_1_l2[0], model_1_linf[0], model_1_knorm[0]),
        (model_1_l1[1], model_1_l2[1], model_1_linf[1], model_1_knorm[1]),
        (model_2_l1[0], model_2_l2[0], model_2_linf[0], model_2_knorm[0]),
        (model_2_l1[1], model_2_l2[1], model_2_linf[1], model_2_knorm[1]),
    ]

    mechanisms = ["l1", "l2", "linf", "Optimal K-norm"]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    y_min = min(min(min(data) for data in data_set) for data_set in data_sets)
    y_max = max(max(max(data) for data in data_set) for data_set in data_sets) + 20

    for i, ax in enumerate(axes.flatten()):
        for data, color, mech in zip(data_sets[i], colors, mechanisms):
            ax.plot(epsilon_values, data, marker="o", linestyle="-", color=color, label=mech)
        ax.set_title(titles[i])
        ax.set_xlabel(r"$\epsilon$")
        ax.set_ylabel("L2 Cost")
        ax.grid(True)
        ax.legend()
        ax.set_ylim([y_min, y_max])

    plt.tight_layout()
    plt.savefig(PLOTS / "Knorm_comparison.png")
    plt.close(fig)


def main() -> None:
    make_gaussian_plot()
    make_knorm_plot()
    print(f"Wrote {PLOTS / 'gaussian_comparison.png'}")
    print(f"Wrote {PLOTS / 'Knorm_comparison.png'}")


if __name__ == "__main__":
    main()
