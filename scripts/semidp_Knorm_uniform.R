# Comparing K-norm mechanisms on Contingency table from uniform probability 
#rm(list=ls())

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
idx <- grep(paste0("^", file_arg), args_file)
if (length(idx) == 0) {
  root_dir <- normalizePath(getwd(), mustWork = FALSE)
} else {
  script_path <- normalizePath(sub(file_arg, "", args_file[idx[1]]), mustWork = FALSE)
  root_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
}
output_dir <- file.path(root_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load required packages
library(MASS)
library(lpSolve)
library(LaplacesDemon)

# Experiment setup
k_values <- c(2,3)
epsilon_values <- c(0.1, 0.5, 1)
num_simulations <- 30
set.seed(123)  # Set global seed for reproducibility
seeds <- sample.int(10000, num_simulations)
n <- 500  # Sample size
SensL1 <- 2
SensL2 <- sqrt(2)
SensLinf <- 1
SensK <- 1
SemiAdj <- 3  


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


# Function to find the basis vectors from a set of vectors
find_basis_vectors <- function(vectors) {
  qr_decomp <- qr(vectors)
  rank <- qr_decomp$rank
  basis <- vectors[, qr_decomp$pivot[1:rank], drop = FALSE]
  return(basis)
}

# Function to check if a vector is a convex combination of a set of vectors
is_convex_combination <- function(v, set_of_vectors) {
  # Convert inputs to numeric and ensure v is a numeric vector
  v <- as.numeric(v)
  set_of_vectors <- as.matrix(set_of_vectors)
  
  # Number of vectors and dimensions
  d <- nrow(set_of_vectors)
  n <- ncol(set_of_vectors)
  
  # Objective function: not minimizing or maximizing anything, so set to zeros
  f.obj <- rep(0, n)
  
  # Constraints: Coefficients must sum to 1 (convex combination)
  # and the combination must equal the vector v
  f.con <- rbind(rep(1, n), set_of_vectors)
  f.dir <- c("=", rep("=", d))
  f.rhs <- c(1, v)
  
  # Solve the linear program
  solution <- lp("min", f.obj, f.con, f.dir, f.rhs)
  
  # Check if a solution was found
  if (solution$status == 0) {
    # Verify if the solution matches the vector v within tolerance
    solution_vector <- set_of_vectors %*% solution$solution
    return(all(abs(solution_vector - v) < 1e-6))
  } else {
    return(FALSE)
  }
}

# Generate K-norm noise
Knorm_noise <- function(eps, SensK, points, max_iterations = 10000) {
  basis_vectors <- find_basis_vectors(points)
  s <- ncol(basis_vectors)
  noise <- NULL
  iteration <- 0
  
  repeat {
    iteration <- iteration + 1
    print(iteration)
    # Generate s:rank iid random variables from U(-1, 1)
    U <- runif(s, min = -1, max = 1)
    # Compute the linear combination of basis vectors
    V <- basis_vectors %*% U
    
    # Check if the noise_candidate is within the convex hull of the points
    if (is_convex_combination(V, points)) {
      noise <- V
      break
    }
    
    # Stop if maximum iterations are reached
    if (iteration >= max_iterations) {
      stop("Maximum iterations reached without finding a valid noise vector.")
    }
    
    # Print a message every 500 iterations
    #if (iteration %% 500 == 0) {
    #  cat("Iteration:", iteration, "- Still searching...\n")
    #}
  }
  
  # Draw a random variable r from Gamma(k, eps/SensK)
  r <- rgamma(1, shape = s, rate = eps / SensK)
  noise <- r * noise
  
  return(noise)
}

#
runsimul <- function(eps, k, SensL1, SensL2, SensLinf, SensK, SemiAdj, n, seed){
  # Set seed for reproducibility
  set.seed(seed)
  
  # Generate a multinomial random variable
  p <- generate_equal_vector(k^2) # For k x k table, we need k^2 probabilities
  table <- rmultinom(1, n, p)
  
  # Generate a projection matrix
  config_matrix <- generate_configurations(k)
  
  # Loss for l1-mechanism
  L1Noise <- rlaplace(k^2, 0, (SemiAdj*SensL1)/eps)
  L1 <- table + L1Noise
  lossL1 <- norm(as.vector(L1 - table), type = "2")
  
  # Loss for l2-mechanism
  L2Noise <- mvrnorm(n = 1, mu = rep(0, k^2), Sigma = diag(k^2))
  r <- rgamma(1, shape = k^2, rate = eps / (SemiAdj*SensL2))
  L2Noise <- r * L2Noise/sqrt(sum(L2Noise^2))
  L2 <- table + L2Noise
  lossL2 <- norm(as.vector(L2-table), type = "2")
  
  # Loss for linf-mechanism 
  LinfNoise <- runif(k^2, min = -1, max = 1)*rgamma(1, shape = k^2, rate = eps / (SemiAdj*SensLinf))
  Linf <- table + LinfNoise
  lossLinf <- norm(as.vector(Linf - table), type = "2")
  
  # Loss for K-norm mechanism
  KnormNoise <- Knorm_noise(eps, SensK, config_matrix,max_iterations = 10000)
  Knorm <- table +KnormNoise
  lossKnorm <- norm(as.vector(Knorm - table), type = "2")
  
  # Return the losses as a named vector
  return(c(LossL1 = lossL1, LossL2 = lossL2, LossLinf = lossLinf, LossKnorm = lossKnorm))
}

# Store results in a data frame
results <- data.frame(k = integer(), epsilon = numeric(), 
                      LossL1 = numeric(), LossL2 = numeric(), 
                      LossLinf = numeric(), LossKnorm = numeric())

# Run the simulations
for (k in k_values) {
  for (eps in epsilon_values) {
    # Initialize vectors to store the losses
    lossL1s <- numeric(num_simulations)
    lossL2s <- numeric(num_simulations)
    lossLinfs <- numeric(num_simulations)
    lossKnorms <- numeric(num_simulations)
    
    for (i in 1:num_simulations) {
      seed <- seeds[i]
      result <- runsimul(eps, k, SensL1, SensL2, SensLinf, SensK, SemiAdj, n, seed)
      lossL1s[i] <- result["LossL1"]
      lossL2s[i] <- result["LossL2"]
      lossLinfs[i] <- result["LossLinf"]
      lossKnorms[i] <- result["LossKnorm"]
    }
    
    # Calculate averages
    avg_lossL1 <- mean(lossL1s)
    avg_lossL2 <- mean(lossL2s)
    avg_lossLinf <- mean(lossLinfs)
    avg_lossKnorm <- mean(lossKnorms)
    
    # Store the results
    results <- rbind(results, data.frame(k = k, epsilon = eps, 
                                         LossL1 = avg_lossL1, LossL2 = avg_lossL2, 
                                         LossLinf = avg_lossLinf, LossKnorm = avg_lossKnorm))
  }
}

# Print the results
print(results)

write.csv(results, file = file.path(output_dir, "Knorm_uniform.csv"))
