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

out_csv <- file.path(results_dir, "grouped_regression_mu_gdp.csv")
out_png <- file.path(plots_dir, "grouped_regression_mu_gdp.png")

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

make_groups <- function(K) {
  grid <- expand.grid(rep(list(c(-1, 1)), K))
  as.matrix(grid)
}

make_design <- function(K) {
  groups <- make_groups(K)
  Z <- cbind(1, groups) / sqrt(K + 1)
  rownames(Z) <- apply(groups, 1, paste, collapse = "")
  Z
}

balanced_counts <- function(N, J) {
  if (N %% J != 0) {
    stop("N must be divisible by J for the balanced-count design.")
  }
  rep(N / J, J)
}

gram_matrix <- function(Z, counts) {
  d <- ncol(Z)
  G <- matrix(0, d, d)
  for (j in seq_len(nrow(Z))) {
    zj <- matrix(Z[j, ], ncol = 1)
    G <- G + counts[j] * (zj %*% t(zj))
  }
  G
}

sensitivities <- function(Z) {
  J <- nrow(Z)
  Lz <- max(sqrt(rowSums(Z^2)))
  delta_b_pair <- 0
  delta_gamma <- 0
  for (j in seq_len(J)) {
    for (k in seq_len(J)) {
      delta_b_pair <- max(delta_b_pair, sqrt(sum((Z[k, ] - Z[j, ])^2)))
      zj <- matrix(Z[j, ], ncol = 1)
      zk <- matrix(Z[k, ], ncol = 1)
      A <- zk %*% t(zk) - zj %*% t(zj)
      delta_gamma <- max(delta_gamma, sqrt(sum(A^2)))
    }
  }
  list(
    Lz = Lz,
    delta_b_semi = 2 * Lz,
    delta_b_1 = max(Lz, delta_b_pair),
    delta_gamma_1 = delta_gamma
  )
}

sample_symmetric_noise <- function(d, sigma) {
  Z <- matrix(0, d, d)
  for (a in seq_len(d)) {
    Z[a, a] <- rnorm(1, mean = 0, sd = sigma)
  }
  if (d >= 2) {
    for (a in seq_len(d - 1)) {
      for (b in (a + 1):d) {
        coef <- rnorm(1, mean = 0, sd = sigma)
        Z[a, b] <- coef / sqrt(2)
        Z[b, a] <- coef / sqrt(2)
      }
    }
  }
  Z
}

project_psd <- function(A) {
  A <- (A + t(A)) / 2
  eig <- eigen(A, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  eig$vectors %*% (vals * t(eig$vectors))
}

coefficient_loss <- function(beta_private, beta_hat) {
  sqrt(sum((beta_private - beta_hat)^2))
}

prediction_loss <- function(Z, beta_private, beta_hat) {
  diff <- as.vector(Z %*% (beta_private - beta_hat))
  sqrt(mean(diff^2))
}

run_one_replicate <- function(K, N, mu, seed, beta_true) {
  set.seed(seed)
  Z <- make_design(K)
  J <- nrow(Z)
  d <- ncol(Z)
  counts <- balanced_counts(N, J)
  sens <- sensitivities(Z)
  Gamma <- gram_matrix(Z, counts)

  prob <- sigmoid(as.vector(Z %*% beta_true[seq_len(d)]))
  sums <- rbinom(J, size = counts, prob = prob)
  b <- as.vector(t(Z) %*% sums)
  beta_hat <- solve(Gamma, b)

  Zb_semi <- rnorm(d, mean = 0, sd = sens$delta_b_semi / mu)
  beta_semi <- solve(Gamma, b + Zb_semi)

  mu_gamma <- mu / (2 * sqrt(2))
  mu_b <- mu / (2 * sqrt(2))
  sigma_gamma <- sens$delta_gamma_1 / mu_gamma
  sigma_b <- sens$delta_b_1 / mu_b
  Gamma_tilde <- Gamma + sample_symmetric_noise(d, sigma_gamma)
  Gamma_tilde_psd <- project_psd(Gamma_tilde)
  b_tilde <- b + rnorm(d, mean = 0, sd = sigma_b)
  beta_dp <- solve(Gamma_tilde_psd, b_tilde)

  data.frame(
    K = K,
    J = J,
    d = d,
    N = N,
    mu = mu,
    seed = seed,
    method = c("Semi-DP: noisy X'Y, exact X'X",
               "Ordinary baseline: noisy X'X and X'Y"),
    coeff_l2 = c(coefficient_loss(beta_semi, beta_hat),
                 coefficient_loss(beta_dp, beta_hat)),
    pred_l2 = c(prediction_loss(Z, beta_semi, beta_hat),
                prediction_loss(Z, beta_dp, beta_hat)),
    coeff_sq = c(sum((beta_semi - beta_hat)^2),
                 sum((beta_dp - beta_hat)^2)),
    pred_sq = c(mean(as.vector(Z %*% (beta_semi - beta_hat))^2),
                mean(as.vector(Z %*% (beta_dp - beta_hat))^2)),
    Lz = sens$Lz,
    delta_b_semi = sens$delta_b_semi,
    delta_b_1 = sens$delta_b_1,
    delta_gamma_1 = sens$delta_gamma_1,
    mu_gamma = c(NA_real_, mu_gamma),
    mu_b = c(NA_real_, mu_b),
    sigma_gamma = c(NA_real_, sigma_gamma),
    sigma_b = c(sens$delta_b_semi / mu, sigma_b),
    gamma_min_eigen = min(eigen(Gamma, symmetric = TRUE, only.values = TRUE)$values),
    gamma_max_eigen = max(eigen(Gamma, symmetric = TRUE, only.values = TRUE)$values),
    stringsAsFactors = FALSE
  )
}

summarize_results <- function(dat) {
  keys <- c("scenario", "x_value", "K", "J", "d", "N", "mu", "method")
  metric_cols <- c("coeff_l2", "pred_l2", "coeff_sq", "pred_sq")
  split_dat <- split(dat, dat[keys], drop = TRUE)
  do.call(rbind, lapply(split_dat, function(df) {
    out <- df[1, keys, drop = FALSE]
    out$n_reps <- nrow(df)
    for (m in metric_cols) {
      vals <- df[[m]]
      out[[paste0(m, "_mean")]] <- mean(vals)
      out[[paste0(m, "_se")]] <- stats::sd(vals) / sqrt(length(vals))
    }
    for (m in c("Lz", "delta_b_semi", "delta_b_1", "delta_gamma_1",
                "mu_gamma", "mu_b", "sigma_gamma", "sigma_b",
                "gamma_min_eigen", "gamma_max_eigen")) {
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
  par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.6, 1.2), oma = c(0, 0, 0, 0))

  cols <- c("Semi-DP: noisy X'Y, exact X'X" = "#1f77b4",
            "Ordinary baseline: noisy X'X and X'Y" = "#ff7f0e")
  pchs <- c("Semi-DP: noisy X'Y, exact X'X" = 16,
            "Ordinary baseline: noisy X'X and X'Y" = 16)
  ltys <- c("Semi-DP: noisy X'Y, exact X'X" = 1,
            "Ordinary baseline: noisy X'X and X'Y" = 2)

  draw_panel <- function(dat, xlab, main, legend_pos, xaxt = "s", axis_labels = NULL) {
    methods <- names(cols)
    ylim <- range(dat$coeff_l2_mean + dat$coeff_l2_se, dat$coeff_l2_mean - dat$coeff_l2_se)
    ylim[1] <- 0
    plot(NA, xlim = range(dat$x_value), ylim = ylim, xlab = xlab,
         ylab = "L2 Cost",
         main = main, xaxt = xaxt, cex.main = 0.9)
    grid(col = "gray70", lty = "solid")
    if (!is.null(axis_labels)) {
      axis(1, at = sort(unique(dat$x_value)), labels = axis_labels)
    }
    for (meth in methods) {
      dd <- dat[dat$method == meth, ]
      dd <- dd[order(dd$x_value), ]
      lines(dd$x_value, dd$coeff_l2_mean, col = cols[meth], lwd = 2, lty = ltys[meth])
      points(dd$x_value, dd$coeff_l2_mean, col = cols[meth], pch = pchs[meth])
      arrows(dd$x_value, dd$coeff_l2_mean - dd$coeff_l2_se,
             dd$x_value, dd$coeff_l2_mean + dd$coeff_l2_se,
             angle = 90, code = 3, length = 0.04, col = cols[meth])
    }
    legend(legend_pos, legend = c("Semi-DP", "Ordinary baseline"),
           col = cols, lty = ltys, pch = pchs, lwd = 2, bty = "n")
  }

  design_dat <- summary_dat[summary_dat$scenario == "design_complexity", ]
  labels <- paste0("K=", sort(unique(design_dat$K)))
  draw_panel(
    design_dat,
    xlab = "Number of binary factors",
    main = "Design complexity trend, N = 512, mu = 1",
    legend_pos = "topleft",
    xaxt = "n",
    axis_labels = labels
  )

  sample_dat <- summary_dat[summary_dat$scenario == "sample_size", ]
  draw_panel(
    sample_dat,
    xlab = "Sample size N",
    main = "Sample-size trend, K = 4, mu = 1",
    legend_pos = "topright"
  )
}

beta_base <- c(0.1, 0.7, -0.5, 0.4, -0.3)
seeds <- seq_len(30)
mu_target <- 1

all_results <- list()

idx <- 1
for (K in 1:4) {
  d <- K + 1
  beta_true <- beta_base[seq_len(d)]
  for (seed in seeds) {
    res <- run_one_replicate(K = K, N = 512, mu = mu_target,
                             seed = 1000 * K + seed, beta_true = beta_true)
    res$scenario <- "design_complexity"
    res$x_value <- K
    all_results[[idx]] <- res
    idx <- idx + 1
  }
}

for (N in c(128, 256, 512, 1024)) {
  K <- 4
  beta_true <- beta_base
  for (seed in seeds) {
    res <- run_one_replicate(K = K, N = N, mu = mu_target,
                             seed = 10000 + N + seed, beta_true = beta_true)
    res$scenario <- "sample_size"
    res$x_value <- N
    all_results[[idx]] <- res
    idx <- idx + 1
  }
}

raw_dat <- do.call(rbind, all_results)
summary_dat <- summarize_results(raw_dat)
summary_dat <- summary_dat[order(summary_dat$scenario, summary_dat$x_value, summary_dat$method), ]
write.csv(summary_dat, out_csv, row.names = FALSE)
plot_results(summary_dat, out_png)

cat("Wrote", out_csv, "\n")
cat("Wrote", out_png, "\n")
