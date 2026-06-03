#!/usr/bin/env Rscript

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
idx <- grep(paste0("^", file_arg), args_file)
if (length(idx) == 0) {
  stop("Run this script with Rscript so that the script path is available.")
}
script_path <- normalizePath(sub(file_arg, "", args_file[idx[1]]), mustWork = FALSE)
root_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

results_dir <- file.path(root_dir, "results")
plots_dir <- file.path(root_dir, "plots")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

out_csv <- file.path(results_dir, "grouped_anova_mu_gdp.csv")
out_png <- file.path(plots_dir, "grouped_anova_mu_gdp.png")

clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

project_simplex <- function(y, total) {
  if (total < 0) {
    stop("Simplex total must be nonnegative.")
  }
  u <- sort(y, decreasing = TRUE)
  cssv <- cumsum(u)
  idx <- seq_along(u)
  active <- which(u - (cssv - total) / idx > 0)
  if (length(active) == 0) {
    rep(total / length(y), length(y))
  } else {
    rho <- max(active)
    theta <- (cssv[rho] - total) / rho
    pmax(y - theta, 0)
  }
}

mean_l2_cost <- function(theta_private, theta_hat) {
  sqrt(mean((theta_private - theta_hat)^2))
}

run_one_replicate <- function(config_name, counts, mu, data_seed, noise_seed) {
  set.seed(data_seed)
  J <- length(counts)
  N <- sum(counts)
  prob <- 0.2 + 0.6 * (seq_len(J) - 1) / (J - 1)

  sums <- rbinom(J, size = counts, prob = prob)
  theta_hat <- sums / counts

  delta_s_semi <- 2
  delta_n_1 <- sqrt(2)
  delta_s_1 <- sqrt(2)

  set.seed(noise_seed)
  s_tilde_semi <- sums + rnorm(J, mean = 0, sd = delta_s_semi / mu)
  theta_semi <- clip01(s_tilde_semi / counts)

  mu_n <- mu / (2 * sqrt(2))
  mu_s <- mu / (2 * sqrt(2))
  n_tilde <- counts + rnorm(J, mean = 0, sd = delta_n_1 / mu_n)
  s_tilde_dp <- sums + rnorm(J, mean = 0, sd = delta_s_1 / mu_s)
  n_plus <- project_simplex(n_tilde, total = N)
  theta_dp <- clip01(s_tilde_dp / pmax(n_plus, 1))

  data.frame(
    config = config_name,
    J = J,
    N = N,
    mu = mu,
    data_seed = data_seed,
    noise_seed = noise_seed,
    method = c("Semi-DP: noisy sums, exact counts",
               "Ordinary baseline: noisy counts and sums"),
    mean_l2 = c(mean_l2_cost(theta_semi, theta_hat),
                mean_l2_cost(theta_dp, theta_hat)),
    mean_sq = c(mean((theta_semi - theta_hat)^2),
                mean((theta_dp - theta_hat)^2)),
    min_count = min(counts),
    max_count = max(counts),
    delta_s_semi = delta_s_semi,
    delta_n_1 = c(NA_real_, delta_n_1),
    delta_s_1 = c(NA_real_, delta_s_1),
    mu_n = c(NA_real_, mu_n),
    mu_s = c(NA_real_, mu_s),
    sigma_n = c(NA_real_, delta_n_1 / mu_n),
    sigma_s = c(delta_s_semi / mu, delta_s_1 / mu_s),
    stringsAsFactors = FALSE
  )
}

summarize_results <- function(dat) {
  keys <- c("config", "J", "N", "mu", "method")
  metric_cols <- c("mean_l2", "mean_sq")
  split_dat <- split(dat, dat[keys], drop = TRUE)
  do.call(rbind, lapply(split_dat, function(df) {
    out <- df[1, keys, drop = FALSE]
    out$n_reps <- nrow(df)
    for (m in metric_cols) {
      vals <- df[[m]]
      out[[paste0(m, "_mean")]] <- mean(vals)
      out[[paste0(m, "_se")]] <- stats::sd(vals) / sqrt(length(vals))
    }
    for (m in c("min_count", "max_count", "delta_s_semi", "delta_n_1",
                "delta_s_1", "mu_n", "mu_s", "sigma_n", "sigma_s")) {
      out[[m]] <- mean(df[[m]], na.rm = TRUE)
    }
    out
  }))
}

plot_results <- function(summary_dat, out_png) {
  png(out_png, width = 1400, height = 600, res = 120)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  par(mfrow = c(1, 2), mar = c(4.6, 5.6, 3.6, 1.2), oma = c(0, 0, 0, 0))

  cols <- c("Semi-DP: noisy sums, exact counts" = "#1f77b4",
            "Ordinary baseline: noisy counts and sums" = "#ff7f0e")
  pchs <- c("Semi-DP: noisy sums, exact counts" = 16,
            "Ordinary baseline: noisy counts and sums" = 16)
  ltys <- c("Semi-DP: noisy sums, exact counts" = 1,
            "Ordinary baseline: noisy counts and sums" = 2)

  draw_panel <- function(config_name, main) {
    dat <- summary_dat[summary_dat$config == config_name, ]
    methods <- names(cols)
    ylim <- range(dat$mean_l2_mean + dat$mean_l2_se,
                  dat$mean_l2_mean - dat$mean_l2_se)
    ylim[1] <- 0
    plot(NA, xlim = range(dat$mu), ylim = ylim, xlab = "mu",
         ylab = "L2 Cost", main = main, cex.main = 0.9)
    grid(col = "gray70", lty = "solid")
    for (meth in methods) {
      dd <- dat[dat$method == meth, ]
      dd <- dd[order(dd$mu), ]
      lines(dd$mu, dd$mean_l2_mean, col = cols[meth], lwd = 2, lty = ltys[meth])
      points(dd$mu, dd$mean_l2_mean, col = cols[meth], pch = pchs[meth])
      arrows(dd$mu, dd$mean_l2_mean - dd$mean_l2_se,
             dd$mu, dd$mean_l2_mean + dd$mean_l2_se,
             angle = 90, code = 3, length = 0.04, col = cols[meth])
    }
    legend("topright", legend = c("Semi-DP", "Ordinary baseline"),
           col = cols, lty = ltys, pch = pchs, lwd = 2, bty = "n")
  }

  draw_panel("balanced", "Balanced counts")
  draw_panel("imbalanced", "Imbalanced counts")
}

counts_list <- list(
  balanced = rep(50, 10),
  imbalanced = c(5, 10, 20, 30, 40, 50, 60, 80, 90, 115)
)
mu_grid <- c(0.25, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0)
seeds <- seq_len(30)

all_results <- list()
idx <- 1
for (config_name in names(counts_list)) {
  config_offset <- if (config_name == "balanced") 100000 else 200000
  for (seed in seeds) {
    data_seed <- config_offset + seed
    for (m_idx in seq_along(mu_grid)) {
      mu <- mu_grid[m_idx]
      noise_seed <- config_offset + 1000 * seed + m_idx
      all_results[[idx]] <- run_one_replicate(
        config_name = config_name,
        counts = counts_list[[config_name]],
        mu = mu,
        data_seed = data_seed,
        noise_seed = noise_seed
      )
      idx <- idx + 1
    }
  }
}

raw_dat <- do.call(rbind, all_results)
summary_dat <- summarize_results(raw_dat)
summary_dat <- summary_dat[order(summary_dat$config, summary_dat$mu,
                                 summary_dat$method), ]
write.csv(summary_dat, out_csv, row.names = FALSE)
plot_results(summary_dat, out_png)

cat("Wrote", out_csv, "\n")
cat("Wrote", out_png, "\n")
