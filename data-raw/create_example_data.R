# Create example dataset for CMP package

set.seed(123)
n <- 500

# Generate predictors
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rbinom(n, 1, 0.5)

# Create correlated error terms
library(MASS)
Sigma_error <- matrix(c(1, 0.5, 0.3, 0.2,
                       0.5, 1, 0.4, 0.1,
                       0.3, 0.4, 1, 0.6,
                       0.2, 0.1, 0.6, 1), nrow = 4)

errors <- mvrnorm(n, mu = rep(0, 4), Sigma = Sigma_error)

# Generate outcomes with correlation through error terms
# Continuous outcome
y_continuous <- 1 + 0.5 * x1 - 0.3 * x2 + 0.8 * x3 + errors[, 1]

# Binary outcome (probit)
y_binary_latent <- -0.2 + 0.7 * x1 + 0.4 * x2 + errors[, 2]
y_binary <- as.numeric(y_binary_latent > 0)

# Ordered outcome (3 categories)
y_ordered_latent <- 0.3 * x1 - 0.5 * x2 + 0.6 * x3 + errors[, 3]
y_ordered <- cut(y_ordered_latent, breaks = c(-Inf, -0.5, 0.5, Inf), labels = 1:3)
y_ordered <- as.numeric(y_ordered)

# Left-censored outcome (Tobit)
y_censored_latent <- 0.8 + 0.4 * x1 + 0.6 * x2 - 0.3 * x3 + errors[, 4]
y_censored <- pmax(0, y_censored_latent)

# Create final dataset
cmp_example <- data.frame(
  id = 1:n,
  x1 = x1,
  x2 = x2,
  x3 = x3,
  y_continuous = y_continuous,
  y_binary = y_binary,
  y_ordered = y_ordered,
  y_censored = y_censored
)

# Save dataset
usethis::use_data(cmp_example, overwrite = TRUE)