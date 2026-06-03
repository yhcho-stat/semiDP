# Gaussian Mechanism on Contingency table from uniform probability 
#rm(list = ls())

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
idx <- grep(paste0("^", file_arg), args_file)
if (length(idx) == 0) {
  root_dir <- normalizePath(getwd(), mustWork = FALSE)
} else {
  script_path <- normalizePath(sub(file_arg, "", args_file[idx[1]]), mustWork = FALSE)
  root_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
}
results_dir <- file.path(root_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Load the required packages
library(pracma)
library(MASS)

# Hyperparameters
mu <- 0.1  # Privacy parameter for the Gaussian mechanism 
n <- 500  # Number of samples in the contingency table
Sens <- sqrt(2)
SensSemi <- 2
SemiAdj <- 3
num_simulations <- 30

# Set the global seed for reproducibility
set.seed(123)
# Generate random seeds for each simulation
seeds <- sample.int(10000, num_simulations)

# Functions to generate projection matrix with different size of contingency table
generate_configurations <- function(k) {
  # Initialize an empty list to store the vectors
  vectors <- list()
  
  # Loop through all pairs of positions for 1s and -1s
  for (i in 1:k) {
    for (j in 1:k) {
      for (m in 1:k) {
        for (n in 1:k) {
          # Ensure the positions are distinct and preserve marginal counts
          if (i != m && j != n) {
            vector <- matrix(0, nrow = k, ncol = k)
            vector[i, j] <- 1
            vector[m, n] <- 1
            vector[i, n] <- -1
            vector[m, j] <- -1
            vectors <- append(vectors, list(as.vector(vector)))
          }
        }
      }
    }
  }
  
  # Convert the list of vectors to a matrix and remove duplicates
  result_matrix <- do.call(cbind, vectors)
  unique_matrix <- unique(t(result_matrix))
  
  return(t(unique_matrix))
}

# Function to find orthonormal columns
find_orthonormal_columns <- function(vectors) {
  vectors <- as.matrix(vectors)
  qr_decomp <- qr(vectors)
  rank <- qr_decomp$rank
  independent_columns <- vectors[, qr_decomp$pivot[1:rank], drop = FALSE]
  Q <- gramSchmidt(independent_columns)$Q
  
  # Ensure Q is orthonormal
  if (any(abs(colSums(Q * Q) - 1) > 1e-10)) {
    stop("Orthonormalization failed: Q columns are not normalized.")
  }
  
  return(Q)
}

# Function to compute the projection matrix
projection_matrix <- function(vectors) {
  Q <- find_orthonormal_columns(vectors)
  P <- Q %*% t(Q)
  return(P)
}

# Function to generate a vector of equal probabilities for a contingency table
generate_equal_vector <- function(k) {
  # Check if k is a positive integer
  if (k <= 0 || k != as.integer(k)) {
    stop("k must be a positive integer")
  }
  
  # Create a vector of length k with each element equal to 1/k
  result_vector <- rep(1/k, k)
  
  return(result_vector)
}

# Function to generate a vector of linearly increasing probabilities for a contingency table
generate_linear_vector <- function(k) {
  # Calculate the total number of cells
  total_cells <- k^2
  
  # Create a vector with linearly increasing values
  values <- 1:total_cells
  
  # Normalize the values to sum to 1
  probabilities <- values / sum(values)
  
  return(probabilities)
}

# Function to generate noise with the covariance structure (SensSemi / mu)^2 * P
generate_noise <- function(P, SensSemi, mu) {
  # Eigenvalue decomposition of P
  eigen_decomp <- eigen(P)
  U <- eigen_decomp$vectors
  Lambda <- diag(eigen_decomp$values)
  
  # Handle small or negative eigenvalues by setting them to zero
  Lambda[Lambda < 1e-10] <- 0
  
  # Compute the noise, considering only the positive part of Lambda
  sqrt_Lambda <- sqrt(Lambda)
  noise <- (SensSemi / mu) * U %*% sqrt_Lambda %*% rnorm(ncol(U))
  
  return(noise)
}

runsimul <- function(mu, k, Sens, SensSemi, SemiAdj, n, seed) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Generate a multinomial random variable
  p <- generate_equal_vector(k^2) # For k x k table, we need k^2 probabilities
  table <- rmultinom(1, n, p)
  
  # Generate a projection matrix
  config_matrix <- generate_configurations(k)
  P <- projection_matrix(config_matrix)
  
  # Generate noise with the correct covariance structure
  noise <- generate_noise(P, SensSemi, mu)
  
  # Compute the 'SemiDP' variable and its loss
  SemiDP <- table + noise
  SemiDP <- pmax(SemiDP, 0)
  lossSemiDP <- norm(as.vector(SemiDP - table), type = "2")
  
  # Compute the 'naive' variable and its loss
  Naive <- table + (SemiAdj * Sens / mu) * mvrnorm(1, numeric(k^2), diag(k^2))
  Naive <- pmax(Naive, 0)
  lossNaive <- norm(as.vector(Naive - table), type = "2")
  
  # Return the losses as a named vector
  return(c(LossSemiDP = lossSemiDP, LossNaive = lossNaive))
}

# Function to run the simulation multiple times and calculate averages
run_multiple_simulations <- function(num_simulations, mu, k, Sens, SensSemi, SemiAdj, n, seeds) {
  results <- lapply(seeds, function(seed) {
    runsimul(mu, k, Sens, SensSemi, SemiAdj, n, seed)
  })
  avg_results <- rowMeans(do.call(cbind, results))
  return(avg_results)
}

# Initialize a list to store the average losses for each k
all_average_losses_uniform <- data.frame()

# Loop over k from 2 to 
for (k in c(2:10)) {
  # Run the simulations and get the average losses
  average_losses_uniform <- run_multiple_simulations(num_simulations, mu, k, Sens, SensSemi, SemiAdj, n, seeds)
  # Store the results in the dataframe
  all_average_losses_uniform <- rbind(all_average_losses_uniform, 
                                      data.frame(k = k, LossSemiDP = average_losses_uniform[1], 
                                                 LossNaive = average_losses_uniform[2]))
}

# Remove row names
rownames(all_average_losses_uniform) <- NULL

# Sort the results by k
all_average_losses_uniform <- all_average_losses_uniform[order(all_average_losses_uniform$k), ]

# Print the results
print(all_average_losses_uniform)
write.csv(all_average_losses_uniform, file = file.path(results_dir, "gaussian_uniform.csv"))
